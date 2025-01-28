import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"
run_out_dir = f"{output_dir}/{{sample_ID}}/{{run_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]

primary_directory = "/home/sak0914/Mtb_Megapipe"


rule create_lineage_helper_files:
    input:
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    params:
        lineage_pos_for_F2 = os.path.join(primary_directory, references_dir, "phylogeny", "Coll2014_positions_all.txt"),
        output_dir = output_dir,
    output:
        bcf_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf",
        bcf_index_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf.csi",
        vcf_lineage_positions = f"{sample_out_dir}/lineage/{{sample_ID}}_lineage_positions.vcf",
    conda:
        f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        # convert the full VCF file to a BCF fileto get only the lineage-defining positions according to the Coll 2014 scheme
        bcftools view {input.vcf_file_gzip} -O b -o {output.bcf_file}

        # index bcf file
        bcftools index {output.bcf_file}

        # create VCF file of just the lineage positions, which will be used by the F2 metric script. Per the documentation, if --regions-file is a tab-delimited file, then it needs two columns (CHROM and POS), and POS is 1-indexed and inclusive
        # THIS IS DIFFERENT BEHAVIOR FROM IF IT WAS A BED FILE OR IF YOU USE BEDTOOLS. IN BOTH OF THOSE CASES, YOU NEED THREE COLUMNS (CHROM, BEG, AND END), AND THEY ARE 0-INDEXED WITH END BEING EXCLUSIVE (I.E. HALF-OPEN)
        bcftools view {output.bcf_file} --regions-file {params.lineage_pos_for_F2} -O v -o {output.vcf_lineage_positions}   
        """


rule lineage_typing:
    input:
        bcf_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf",
        bcf_index_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf.csi",
        vcf_lineage_positions = f"{sample_out_dir}/lineage/{{sample_ID}}_lineage_positions.vcf",
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    params:
        lineage_SNP_info = os.path.join(primary_directory, references_dir, "phylogeny", "Coll2014_SNPs_all.csv"),
        F2_metric_script = os.path.join(primary_directory, scripts_dir, "calculate_F2_metric.py"),
        output_dir = output_dir,        
    output:
        F2_metric_output = f"{sample_out_dir}/lineage/F2_Coll2014.txt",
        minor_allele_fractions_output = f"{sample_out_dir}/lineage/minor_allele_fractions.csv",
        fast_lineage_caller_output = f"{sample_out_dir}/lineage/fast_lineage_caller_output.txt",
    shell:
        """
        python3 -u {params.F2_metric_script} -i {params.output_dir}/{wildcards.sample_ID} -o {output.F2_metric_output} -O {output.minor_allele_fractions_output} --lineage-file {params.lineage_SNP_info}

        rm {input.bcf_file} {input.bcf_index_file} {input.vcf_lineage_positions}

        fast-lineage-caller {input.vcf_file_variants_only} --pass --out {output.fast_lineage_caller_output}
        """


rule combine_codon_variants:
    input:
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    output:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.vcf",
    params:
        combine_codon_variants_script = os.path.join(primary_directory, scripts_dir, "combine_codon_variants.py"),
    shell:
        """
        # this creates _variants_combined_codons.vcf in the same directory as the input VCF. You can also specify the -o flag if you want to write the output file somewhere else
        python3 -u {params.combine_codon_variants_script} -i {input.vcf_file_variants_only}        
        """



rule annotate_variants_snpEff:
    input:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_combinedCodons.vcf",
    output:
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf",
    conda:
        f"{primary_directory}/envs/variant_annotation.yaml"
    params:
        snpEff_db = config['snpEff_db'],
    shell:
        """
        chrom_name=$(zgrep -v '^#' {input.vcf_file_variants_combinedCodons} | awk 'NR==1 {{print $1}}')
    
        if [ "$chrom_name" = "Chromosome" ]; then
            snpEff_db="Mycobacterium_tuberculosis_h37rv"
        else
            snpEff_db="Mycobacterium_tuberculosis_gca_000195955"
        fi

        snpEff eff $snpEff_db -noStats -no-downstream -no-upstream -lof {input.vcf_file_variants_combinedCodons} > {output.vcf_file_variants_combinedCodons_annot}

        rm {input.vcf_file_variants_combinedCodons}
        """



rule create_WHO_catalog_variants_TSV:
    input:
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf",
    output:
        vcf_file_bgzip = temp(f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf.bgz"),
        vcf_file_bgzip_tbi = temp(f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_combinedCodons.eff.vcf.bgz.tbi"),
        variants_file_tsv = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants.tsv",
    params:
        # WHO_catalog_regions_BED_file = os.path.join(primary_directory, references_dir, "WHO_catalog_resistance", "regions.bed"),
        WHO_catalog_regions_BED_file = os.path.join(primary_directory, references_dir, "WHO_catalog_resistance", "regions_refseq.bed"),
    conda:
        f"{primary_directory}/envs/variant_annotation.yaml"
    shell:
        """
        chrom_name=$(zgrep -v '^#' {input.vcf_file_variants_combinedCodons_annot} | awk 'NR==1 {{print $1}}')
    
        if [ "$chrom_name" = "Chromosome" ]; then
            bed_file="/home/sak0914/Mtb_Megapipe/references/WHO_catalog_resistance/regions.bed"
        else
            bed_file="/home/sak0914/Mtb_Megapipe/references/WHO_catalog_resistance/regions_refseq.bed"
        fi

        # need to bgzip the VCF file to use bcftools view with the region argument. NEED TO PUT "" AROUND FILE NAME TO PROPERLY CONSIDER SPECIAL CHARACTERS IN FILENAME
        bgzip -c {input.vcf_file_variants_combinedCodons_annot} > {output.vcf_file_bgzip}
    
        # tabix the bgzipped file, which will create fName.bgz.tbi
        tabix -p vcf -f {output.vcf_file_bgzip}

        bcftools view -R $bed_file {output.vcf_file_bgzip} | SnpSift extractFields '-' POS REF ALT FILTER QUAL IMPRECISE AF DP BQ MQ IC DC ANN -e "" > {output.variants_file_tsv}
        """


rule get_WHO_catalog_resistance_predictions:
    input:
        variants_file_tsv = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants.tsv",
    output:
        variants_file_tsv_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_annot.tsv",
        predictions_fName = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_pred_AF_thresh_75.csv",
        predictions_fName_lowAF = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_pred_AF_thresh_25.csv",
    params:
        output_file_basename = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_pred",
        process_variants_WHO_catalog_script = os.path.join(primary_directory, scripts_dir, "process_variants_for_WHO_catalog.py"),
        get_WHO_resistance_predictions_script = os.path.join(primary_directory, scripts_dir, "WHO_catalog_resistance_pred.py"),
    shell:
        """
        python3 -u {params.process_variants_WHO_catalog_script} -i {input.variants_file_tsv}
        rm {input.variants_file_tsv}

        # get resistance predictions -- any Group 1 or 2 variant that passes QC leads to a prediction of R for a given drug. If not, predicted S
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o {params.output_file_basename}
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o {params.output_file_basename} --AF-thresh 0.25
        """