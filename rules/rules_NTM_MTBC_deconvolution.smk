import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"
run_in_dir = f"{input_dir}/{{sample_ID}}/{{run_ID}}"
run_out_dir = f"{output_dir}/{{sample_ID}}/{{run_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]

primary_directory = "/home/sak0914/Mtb_Megapipe"



rule kraken_classification_NTM_MTBC_database:
    input:
        fastq1_trimmed = f"{run_in_dir}/fastp/{{run_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{run_in_dir}/fastp/{{run_ID}}.R2.trimmed.fastq.gz",
    output:
        kraken_report = f"{run_out_dir}/kraken/NTM_MTBC_report.txt",
        kraken_classifications = f"{run_out_dir}/kraken/NTM_MTBC_classifications.txt",
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    params:
        confidence = config['kraken2_confidence'],
        minimum_hit_groups = config['kraken2_minimum_hit_groups'],
        kraken_db = "/n/scratch/users/s/sak0914/mycobacteriaceae",
        output_dir = output_dir,
    shell:
        """
        kraken2 --db {params.kraken_db} --threads 8 \
                --confidence {params.confidence} \
                --minimum-hit-groups {params.minimum_hit_groups} \
                --paired {input.fastq1_trimmed} {input.fastq2_trimmed} \
                --gzip-compressed \
                --report {output.kraken_report} > {output.kraken_classifications}
        
        # rm {input.fastq1_trimmed} {input.fastq2_trimmed}
        """



rule separate_kraken_read_names:
    input:
        kraken_report = f"{run_out_dir}/kraken/NTM_MTBC_report.txt",
        kraken_classifications = f"{run_out_dir}/kraken/NTM_MTBC_classifications.txt",
    output:
        MTBC_read_names = f"{run_out_dir}/kraken/MTBC_read_names.txt",
        BLAST_read_names = f"{run_out_dir}/kraken/BLAST_read_names.txt",
    params:
        separate_kraken_reads_script = os.path.join(primary_directory, "scripts", "separate_kraken_reads.py"),
        kraken_db = "/n/scratch/users/s/sak0914/mycobacteriaceae",
        MTBC_taxid = config['MTBC_taxid'],
    shell:
        """
        python3 -u {params.separate_kraken_reads_script} \
                -i {input.kraken_classifications} \
                -d {params.kraken_db} \
                -t {params.MTBC_taxid} \
        """



rule create_FASTA_for_BLAST:
    input:
        fastq1_trimmed = f"{run_in_dir}/fastp/{{run_ID}}.R1.trimmed.fastq.gz",
        BLAST_read_names = f"{run_out_dir}/kraken/BLAST_read_names.txt",
    output:
        BLAST_FASTA_file = f"{run_out_dir}/BLAST/reads_for_BLAST.fasta"
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    shell:
        """
        # sufficient to just do it for 1 of the read files?
        seqtk subseq {input.fastq1_trimmed} {input.BLAST_read_names} | seqtk seq -a  > {output.BLAST_FASTA_file}
        """



rule BLAST_ambiguous_reads:
    input:
        BLAST_FASTA_file = f"{run_out_dir}/BLAST/reads_for_BLAST.fasta",
    output:
        BLAST_hits_file = f"{run_out_dir}/BLAST/hits.tsv",
    params:
        blast_db = "/n/scratch/users/s/sak0914/mycobacteriaceae_blast_DB/mycobacteriaceae_blast_DB",
    conda:
        # f"{primary_directory}/envs/read_processing_aln.yaml",
        "/home/sak0914/anaconda3/envs/bioinformatics"
    shell:
        """
        # -outfmt 6 means tabular format
        # use megablast, which is good for highly similar sequences (yes because the reads are likely NTM and so will match with high identity to the DB). It's also faster
        blastn -num_threads 8 \
               -max_target_seqs 10 \
               -task megablast \
               -query {input.BLAST_FASTA_file} \
               -db {params.blast_db} \
               -out {output.BLAST_hits_file} \
               -outfmt 6
        """



rule separate_BLAST_reads:
    input:
        BLAST_hits_file = f"{run_out_dir}/BLAST/hits.tsv",
    output:
        BLAST_MTBC_read_names = f"{run_out_dir}/BLAST/MTBC_read_names.txt",        
    params:
        separate_BLAST_reads_script = os.path.join(primary_directory, "scripts", "separate_BLAST_reads.py"),
        kraken_db = "/n/scratch/users/s/sak0914/mycobacteriaceae",
        MTBC_taxid = config['MTBC_taxid'],
    shell:
        """
        python3 -u {params.separate_BLAST_reads_script} \
                -i {input.BLAST_hits_file} \
                -d {params.kraken_db} \
                -t {params.MTBC_taxid} \
        """



rule extract_MTBC_reads:
    input:
        fastq1_trimmed = f"{run_in_dir}/fastp/{{run_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{run_in_dir}/fastp/{{run_ID}}.R2.trimmed.fastq.gz",
        kraken_MTBC_read_names = f"{run_out_dir}/kraken/MTBC_read_names.txt",
        BLAST_MTBC_read_names = f"{run_out_dir}/BLAST/MTBC_read_names.txt",
    output:
        fastq1_trimmed_classified = f"{run_out_dir}/bam/{{run_ID}}.R1.MTBC.fastq",
        fastq2_trimmed_classified = f"{run_out_dir}/bam/{{run_ID}}.R2.MTBC.fastq",
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    shell:
        """
        # concatenate filtered reads from both kraken and BLAST into a single file
        cat <(seqtk subseq {input.fastq1_trimmed} {input.kraken_MTBC_read_names}) <(seqtk subseq {input.fastq1_trimmed} {input.BLAST_MTBC_read_names}) > {output.fastq1_trimmed_classified}
        cat <(seqtk subseq {input.fastq2_trimmed} {input.kraken_MTBC_read_names}) <(seqtk subseq {input.fastq2_trimmed} {input.BLAST_MTBC_read_names}) > {output.fastq2_trimmed_classified}
        """




rule align_reads_mark_duplicates:
    input:
        fastq1_trimmed_classified = f"{run_out_dir}/bam/{{run_ID}}.R1.MTBC.fastq",
        fastq2_trimmed_classified = f"{run_out_dir}/bam/{{run_ID}}.R2.MTBC.fastq",
    output:
        sam_file = temp(f"{run_out_dir}/bam/{{run_ID}}.sam"),
        bam_file = temp(f"{run_out_dir}/bam/{{run_ID}}.bam"),
        bam_index_file = temp(f"{run_out_dir}/bam/{{run_ID}}.bam.bai"),
        bam_file_dedup = f"{run_out_dir}/bam/{{run_ID}}.dedup.bam",
        bam_file_dedup_metrics = f"{run_out_dir}/bam/{{run_ID}}.dedup.bam.metrics",
        bam_index_file_dedup = f"{run_out_dir}/bam/{{run_ID}}.dedup.bam.bai",
    params:
        output_dir = output_dir,
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    shell:
        """
        # index reference genome (which is required before aligning reads)
        bwa-mem2 index {params.ref_genome}

        # align reads to the reference genome sequence. The RG name specifies the read group name, which is necessary if you are merging multiple WGS runs into a single BAM file
        bwa-mem2 mem -M -R "@RG\\tID:{wildcards.run_ID}\\tSM:{wildcards.run_ID}" -t 8 {params.ref_genome} {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} > {output.sam_file}

        # sort alignment and convert to bam file
        samtools view -b {output.sam_file} | samtools sort > {output.bam_file}

        # index alignment, which creates a .bai index file
        samtools index {output.bam_file}

        # -Xmx6g specifies to allocate 6 GB
        picard -Xmx30g MarkDuplicates I={output.bam_file} O={output.bam_file_dedup} REMOVE_DUPLICATES=true M={output.bam_file_dedup_metrics} ASSUME_SORT_ORDER=coordinate READ_NAME_REGEX='(?:.*.)?([0-9]+)[^.]*.([0-9]+)[^.]*.([0-9]+)[^.]*$'

        # index the deduplicated alignment with samtools, which will create a dedup_bam_file.bai file
        samtools index {output.bam_file_dedup}
        """


rule get_BAM_file_depths:
    input:
        bam_file_dedup = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.dedup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]],
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_out_dir = sample_out_dir,
    output:
        depth_file = temp(f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv"),
        depth_file_gzip = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv.gz",
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    shell:
        """
        # Write all .bam files to a text file

        # get all runs associated with this sample_ID and compute depth
        # -a computes depth at all positions, not just those with non-zero depth
        # -Q is for minimum mapping quality: use 1, so that multiply mapped reads aren't counted. These have mapping quality of 0
        samtools depth -a -Q 1 {input.bam_file_dedup} > {output.depth_file}

        # get the length of the reference genome
        genome_length=$(tail -n +2 {params.ref_genome} | tr -d '\n' | wc -c) # remove first line (FASTA header) and newline characters, then count characters to get ref genome length

        # when there are multiple bam files, each one is its own column in the depth file.
        num_sites_H37Rv=$(wc -l {output.depth_file} | awk '{{print $1}}')
    
        if [ ! "$num_sites_H37Rv" -eq "$genome_length" ]; then
            echo "Check that all $genome_length sites in the H37Rv reference genome are in {output.depth_file}, which currently has $num_sites_H37Rv sites"
            exit 1
        fi

        gzip -c {output.depth_file} > {output.depth_file_gzip}
        """


rule get_BAMs_passing_QC_thresholds:
    input:
        depth_file_gzip = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv.gz", # contains depths for all BAM files for all WGS runs
        bam_file_dedup = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.dedup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]],
    output:
        pass_BAMs_file = f"{sample_out_dir}/bam/pass_BAMs.txt",
    params:
        sample_out_dir = sample_out_dir,
        BAM_depth_QC_script = os.path.join(primary_directory, scripts_dir, "BAM_depth_QC.py"),
        median_depth = config["median_depth"],
        min_cov = config["min_cov"],
        genome_cov_prop = config["genome_cov_prop"],
    shell:
        """
        # run the script to determine which runs pass the BAM depth criteria
        python3 -u {params.BAM_depth_QC_script} -i {input.depth_file_gzip} -b {input.bam_file_dedup} --median-depth {params.median_depth} --min-cov {params.min_cov} --genome-cov-prop {params.genome_cov_prop}
        """


rule merge_BAMs:
    input:
        pass_BAMs_file = f"{sample_out_dir}/bam/pass_BAMs.txt",
    output:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
        merged_bam_index_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam.bai",
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_out_dir = sample_out_dir,
        median_depth = config["median_depth"],
        min_cov = config["min_cov"],
        genome_cov_prop = config["genome_cov_prop"],
    shell:
        """
        num_runs_passed=$(wc -l {input.pass_BAMs_file} | awk '{{print $1}}')

        # stop processing samples that don't pass the BAM coverage requirements
        if [ $num_runs_passed -eq 0 ]; then
            echo "No BAM files for {wildcards.sample_ID} passed the minimum coverage requirements. Halting pipeline for this sample"
            exit
            
        else 
            echo "$num_runs_passed WGS runs for {wildcards.sample_ID} have median depth ≥ {params.median_depth} and at least {params.genome_cov_prop} of sites with {params.min_cov}x coverage"

            # merge them using samtools. works because the original bam files were sorted prior to running picard and dropping duplicates (after which they remain sorted)
            samtools merge -b {input.pass_BAMs_file} -o {output.merged_bam_file}

            if [ $num_runs_passed -eq 1 ]; then

                # if only one BAM file passed, delete the original BAM file to reduce disk space usage because it's a duplicate of the merged BAM file
                for file_path in $(cat {input.pass_BAMs_file}); do
                    rm "$file_path" "$file_path.bai"
                done

            fi

            # index the merged BAM file for variant calling
            samtools index {output.merged_bam_file}

        fi
        """


rule pilon_variant_calling:
    input:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    output:
        vcf_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.vcf"),
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
        fasta_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.fasta"),        
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_pilon_dir = f"{sample_out_dir}/pilon",
    conda:
        f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        pilon -Xmx30g --minmq 1 --genome {params.ref_genome} --bam {input.merged_bam_file} --output {wildcards.sample_ID} --outdir {params.sample_pilon_dir} --variant
            
        # then gzip the full VCF file and delete the unzipped version. Also delete the FASTA file because it's not needed
        gzip -c {output.vcf_file} > {output.vcf_file_gzip}

        # save the variants only (non-REF calls) to another VCF file. Left-align indels and deduplicate variants with the same POS, REF, and ALT.
        # this affects those cases where the position of the indel is ambiguous
        # however, because of the shifting positions, the position of the indel can change, so need to sort it again
        bcftools norm --rm-dup none --fasta-ref {params.ref_genome} {output.vcf_file_gzip} | bcftools sort | bcftools view --types snps,indels,mnps,other > {output.vcf_file_variants_only}
        """



rule freebayes_WHO_regions:
    input:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    output:
        vcf_file = temp(f"{sample_out_dir}/freebayes/{{sample_ID}}.vcf"),
        vcf_file_gzip = f"{sample_out_dir}/freebayes/{{sample_ID}}_full.WHO_regions.vcf.gz",
        vcf_file_variants_only = f"{sample_out_dir}/freebayes/{{sample_ID}}.variants.WHO_regions.vcf",
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        WHO_regions = os.path.join(primary_directory, references_dir, "WHO_catalog_resistance", "regions.bed"),
    conda:
        f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        # -p is ploidy
        # freebayes says it automatically does left-alignment of indels, but there was an issue with that in the WHO catalog, so do it as well
        # so left-align indels and drop duplicate records
        freebayes -f {params.ref_genome} -p 1 \
                --min-mapping-quality 30 \
                --min-base-quality 30 \
                --region {params.WHO_regions} \
                -b {input.merged_bam_file} \
                -v {output.vcf_file}

        # then gzip the full VCF file
        gzip -c {output.vcf_file} > {output.vcf_file_gzip}

        # save the variants only (non-REF calls) to another VCF file. Left-align indels and deduplicate variants with the same POS, REF, and ALT.
        # this affects those cases where the position of the indel is ambiguous
        # however, because of the shifting positions, the position of the indel can change, so need to sort it again
        bcftools norm --rm-dup none --fasta-ref {params.ref_genome} {output.vcf_file_gzip} | bcftools sort | bcftools view --types snps,indels,mnps,other > {output.vcf_file_variants_only}
        """



rule freebayes_variant_calling:
    input:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    output:
        vcf_file = temp(f"{sample_out_dir}/freebayes/{{sample_ID}}.vcf"),
        vcf_file_gzip = f"{sample_out_dir}/freebayes/{{sample_ID}}_full.vcf.gz",
        vcf_file_variants_only = f"{sample_out_dir}/freebayes/{{sample_ID}}_variants.vcf",
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
    conda:
        f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        # -p is ploidy
        # freebayes says it automatically does left-alignment of indels, but there was an issue with that in the WHO catalog, so do it as well
        # so left-align indels and drop duplicate records
        # leave --min-alternate-count at the default of 2
        # the minimum AF we're going down to is 1%, so set --min-alternate-fraction to 0.01
        freebayes -f {params.ref_genome} -p 1 --min-mapping-quality 30 --min-base-quality 30 --min-alternate-fraction 0.01 -b {input.merged_bam_file} -v {output.vcf_file}

        # left-align and deduplicate variants with the same POS, REF, and ALT in the full VCF file, then gzip
        bcftools norm --rm-dup none --fasta-ref {params.ref_genome} {output.vcf_file} | bcftools sort | gzip -c > {output.vcf_file_gzip}

        # save the variants only (non-REF calls) to another VCF file. Split multi-allelic sites for easier parsing of the variants
        # however, because of the shifting positions due to left-aligning indels above, the position of the indel can change, so need to sort it again
        bcftools view --types snps,indels,mnps,other {output.vcf_file_gzip} | bcftools norm --multiallelics -any '-' | bcftools sort > {output.vcf_file_variants_only}
        """


rule annotate_extract_freebayes_variants:
    input:
        vcf_file_variants_only = f"{sample_out_dir}/freebayes/{{sample_ID}}_variants.vcf",
    output:
        vcf_file_variants_only_annot = f"{sample_out_dir}/freebayes/{{sample_ID}}_variants.eff.vcf",
        variants_file_tsv = f"{sample_out_dir}/freebayes/{{sample_ID}}_variants.tsv",
        field_names = temp(f"{sample_out_dir}/freebayes/field_names.txt"),
    conda:
        f"{primary_directory}/envs/variant_annotation.yaml"
    params:
        snpEff_db = config['snpEff_db'],
    shell:
        """
        chrom_name=$(zgrep -v '^#' {input.vcf_file_variants_only} | awk 'NR==1 {{print $1}}')
    
        if [ "$chrom_name" = "Chromosome" ]; then
            snpEff_db="Mycobacterium_tuberculosis_h37rv"
        else
            snpEff_db="Mycobacterium_tuberculosis_gca_000195955"
        fi

        # run snpEff annotation
        snpEff eff $snpEff_db -noStats -no-downstream -no-upstream -lof {input.vcf_file_variants_only} > {output.vcf_file_variants_only_annot}

        # get all field names because not sure which ones we will need for low AF variant detection
        echo -e "POS\nREF\nALT\nQUAL\nFILTER\nANN[*].GENE\nANN[*].HGVS_C\nANN[*].HGVS_P\nANN[*].EFFECT" > {output.field_names}
        grep "^##INFO=<ID=" {output.vcf_file_variants_only_annot} | cut -d'=' -f3 | cut -d',' -f1 >> {output.field_names}
        grep "^##FORMAT=<ID=" {output.vcf_file_variants_only_annot} | cut -d'=' -f3 | cut -d',' -f1 >> {output.field_names}
        
        # extract all variants to a TSV file. Keep all fields 
        SnpSift extractFields {output.vcf_file_variants_only_annot} $(paste -sd " " {output.field_names}) -e "" > {output.variants_file_tsv}
        """



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
        minor_allele_fractions_output = temp(f"{sample_out_dir}/lineage/minor_allele_fractions.csv"),
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


rule get_SNPs_for_phylogenetic_tree:
    input:
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    output:
        vcf_SNP_sites = temp(f"{sample_out_dir}/lineage/SNP_sites.tsv"),
        vcf_SNP_sites_gzip = f"{sample_out_dir}/lineage/SNP_sites.tsv.gz",
    conda:
        f"{primary_directory}/envs/variant_annotation.yaml"
    params:
        exclude_regions_BED_file = os.path.join(primary_directory, references_dir, "phylogeny", "exclude_regions.bed"),
    shell:
        """
        # exclude structural variants (SVTYPE)
        gunzip -c {input.vcf_file_gzip} | bcftools filter -e "SVTYPE == 'INS' | SVTYPE == 'DEL'" | bedtools subtract -a '-' -b {params.exclude_regions_BED_file} | SnpSift extractFields '-' POS REF ALT FILTER QUAL AF DP BQ MQ -e "" > {output.vcf_SNP_sites}

        gzip -c {output.vcf_SNP_sites} > {output.vcf_SNP_sites_gzip}
        """