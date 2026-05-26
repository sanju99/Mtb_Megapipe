import os, glob
import numpy as np
import pandas as pd


# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]

conda_directory = config['conda_dir']
primary_directory = "/home/sak0914/MtbLongitudinalDiversity/lowAF_variant_calling"


rule trim_adapters:
    input:
        fastq1 = lambda wildcards: R1_dict[wildcards.sample_ID],
        fastq2 = lambda wildcards: R2_dict[wildcards.sample_ID],
    output:
        fastq1_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}.R2.trimmed.fastq.gz",
        fastp_html = f"{sample_out_dir}/fastp/fastp.html",
        fastp_json = f"{sample_out_dir}/fastp/fastp.json"
    conda:
        f"{conda_directory}/envs/read_processing_aln.yaml"
    params:
        min_read_length = config["min_read_length"]
    threads:
        8
    shell:
        """
        fastp -i {input.fastq1} -I {input.fastq2} \
              -o {output.fastq1_trimmed} -O {output.fastq2_trimmed} \
              -h {output.fastp_html} \
              -j {output.fastp_json} \
              --length_required {params.min_read_length} \
              --dedup \
              --thread {threads}
        """


rule kraken_classification:
    input:
        fastq1_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}.R2.trimmed.fastq.gz",
    output:
        kraken_report = f"{sample_out_dir}/kraken/kraken_report_standard_DB.txt",
        kraken_classifications = f"{sample_out_dir}/kraken/kraken_classifications_standard_DB",
    conda:
        f"{conda_directory}/envs/read_processing_aln.yaml"
    params:
        kraken_db = config['kraken_db'],
        output_dir = output_dir,
    threads:
        8
    shell:
        """
        # --confidence is the minimum fraction of k-mers in a read that must match a given taxon for that read to be assigned to that taxon
        kraken2 --db {params.kraken_db} \
                --threads {threads} \
                --confidence 0 \
                --paired {input.fastq1_trimmed} {input.fastq2_trimmed} \
                --gzip-compressed \
                --report {output.kraken_report} \
                --output {output.kraken_classifications} \
                --memory-mapping
        """
        
        
        
rule extract_kraken_read_names:
    input:
        kraken_classifications = f"{sample_out_dir}/kraken/kraken_classifications_standard_DB",
    output:  
        kraken_classifications_gzipped = f"{sample_out_dir}/kraken/kraken_classifications_standard_DB.csv.gz", # gets gzipped by the python script. Did this to add headers
        keep_read_names = f"{sample_out_dir}/kraken/keep_read_names.txt"
    params:
        kraken_db = config['kraken_db'],
        extract_kraken_reads_script = os.path.join(primary_directory, scripts_dir, "extract_kraken_read_names.py"),
        taxid = config['taxid'],
    shell:
        """
        python3 -u {params.extract_kraken_reads_script} \
                -t {params.taxid} \
                -d {params.kraken_db} \
                -i {input.kraken_classifications} \
                -o {output.keep_read_names} \
                --include-children \
                --include-parents
                
        rm {input.kraken_classifications}
        """



rule extract_kraken_reads:
    input:
        fastq1_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{sample_out_dir}/fastp/{{sample_ID}}.R2.trimmed.fastq.gz",
        keep_read_names = f"{sample_out_dir}/kraken/keep_read_names.txt"
    output:
        fastq1_trimmed_classified = f"{sample_out_dir}/kraken/{{sample_ID}}.R1.kraken.filtered.fastq.gz",
        fastq2_trimmed_classified = f"{sample_out_dir}/kraken/{{sample_ID}}.R2.kraken.filtered.fastq.gz",    
    conda:
        f"{conda_directory}/envs/read_processing_aln.yaml"
    shell:
        """
        # seqtk will write outputs to unzipped files, even if the input was compressed
        seqtk subseq {input.fastq1_trimmed} {input.keep_read_names} | gzip -c > {output.fastq1_trimmed_classified} 
        seqtk subseq {input.fastq2_trimmed} {input.keep_read_names} | gzip -c > {output.fastq2_trimmed_classified} 
        
        rm {input.fastq1_trimmed} {input.fastq2_trimmed}
        
        gzip {input.keep_read_names}
        """
        
        
        
rule fastlin_typing:
    input:
        fastq1_trimmed_classified = f"{sample_out_dir}/kraken/{{sample_ID}}.R1.kraken.filtered.fastq.gz",
        fastq2_trimmed_classified = f"{sample_out_dir}/kraken/{{sample_ID}}.R2.kraken.filtered.fastq.gz",
    output:
        fastq1_trimmed_classified_gzipped = temp(f"{sample_out_dir}/fastlin/{{sample_ID}}_1.fastq.gz"),
        fastq2_trimmed_classified_gzipped = temp(f"{sample_out_dir}/fastlin/{{sample_ID}}_2.fastq.gz"),
        fastlin_dir = directory(f"{sample_out_dir}/fastlin"),
        fastlin_output = f"{sample_out_dir}/fastlin/output.txt",
        fastlin_separate_output = f"{sample_out_dir}/fastlin/sep_output.txt",
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    params:
        kraken_dir = f"{sample_out_dir}/kraken",
        fastlin_barcodes = os.path.join(primary_directory, references_dir, "phylogeny", "MTBC_barcodes.tsv"),
    shell:
        """
        # there's one TRUST sample with coverage > 900
        fastlin -d {params.kraken_dir} -b {params.fastlin_barcodes} -o {output.fastlin_separate_output} -x 1000
        
        # have to do this because fastlin won't recognize that 2 FASTQ files are paired from the same sample unless they have the suffixes _1.fastq.gz and _2.fastq.gz
        cp {input.fastq1_trimmed_classified} {output.fastq1_trimmed_classified_gzipped}
        cp {input.fastq2_trimmed_classified} {output.fastq2_trimmed_classified_gzipped}
        
        # there's one TRUST sample with coverage > 900
        fastlin -d {output.fastlin_dir} -b {params.fastlin_barcodes} -o {output.fastlin_output} -x 1000
        """



rule align_reads_mark_duplicates:
    input:
        fastlin_output = f"{sample_out_dir}/fastlin/output.txt", # so that it's run first because we delete the FASTQ files at the end of this rule
        fastq1_trimmed_classified = f"{sample_out_dir}/kraken/{{sample_ID}}.R1.kraken.filtered.fastq.gz",
        fastq2_trimmed_classified = f"{sample_out_dir}/kraken/{{sample_ID}}.R2.kraken.filtered.fastq.gz",  
    output:
        sam_file = temp(f"{sample_out_dir}/bam/{{sample_ID}}.sam"),
        bam_file = temp(f"{sample_out_dir}/bam/{{sample_ID}}.bam"),
        bam_index_file = temp(f"{sample_out_dir}/bam/{{sample_ID}}.bam.bai"),
        bam_file_dedup = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
        bam_file_dedup_metrics = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam.metrics",
        bam_index_file_dedup = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam.bai",
    params:
        output_dir = output_dir,
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        bwa_mem_seed_length = config['bwa_mem_seed_length']
    conda:
        f"{conda_directory}/envs/read_processing_aln_bwa.yaml"
    threads:
        8
    shell:
        """
        bwa mem -M -R "@RG\\tID:{wildcards.sample_ID}\\tSM:{wildcards.sample_ID}" \
                    -k {params.bwa_mem_seed_length} \
                    -t {threads} \
                    {params.ref_genome} \
                    {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} \
                    > {output.sam_file}

        # sort alignment and convert to bam file
        samtools view -b {output.sam_file} | samtools sort > {output.bam_file}

        # index alignment, which creates a .bai index file
        samtools index {output.bam_file}

        # -Xmx6g specifies to allocate 6 GB
        picard -Xmx10g MarkDuplicates I={output.bam_file} O={output.bam_file_dedup} REMOVE_DUPLICATES=true M={output.bam_file_dedup_metrics} ASSUME_SORT_ORDER=coordinate READ_NAME_REGEX='(?:.*.)?([0-9]+)[^.]*.([0-9]+)[^.]*.([0-9]+)[^.]*$'

        # index the deduplicated alignment with samtools, which will create a dedup_bam_file.bai file
        samtools index {output.bam_file_dedup}
        
        rm {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified}
        """



rule get_BAM_file_depths:
    input:
        bam_file_dedup = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_out_dir = sample_out_dir,
    output:
        depth_file = temp(f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv"),
        depth_file_gzip = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv.gz",
    conda:
        f"{conda_directory}/envs/read_processing_aln.yaml"
    shell:
        """
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

        
        
rule pilon_variant_calling:
    input:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    output:
        vcf_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.vcf"),
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
        fasta_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.fasta"),        
        vcf_SNPs = f"{sample_out_dir}/pilon/{{sample_ID}}_SNPs.vcf",
        vcf_indels = f"{sample_out_dir}/pilon/{{sample_ID}}_indels.vcf",
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_pilon_dir = f"{sample_out_dir}/pilon",
    conda:
        f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        pilon -Xmx10g --minmq 1 --genome {params.ref_genome} --bam {input.merged_bam_file} --output {wildcards.sample_ID} --outdir {params.sample_pilon_dir} --variant
            
        # left-align indels and split multi-allelics, then gzip the full VCF file 
        # this affects those cases where the position of the indel is ambiguous
        bcftools norm --multiallelics -both --fasta-ref {params.ref_genome} {output.vcf_file} | bcftools sort | gzip -c > {output.vcf_file_gzip}
                        
        # save a VCF file of the SNP variants. Anything with BC >= 5 for at least 2 alleles       
        bcftools filter -i '(strlen(REF) == strlen(ALT)) & ((REF=="A" && ((BC[1]>=5) || (BC[2]>=5) || (BC[3]>=5))) || (REF=="C" && ((BC[0]>=5) || (BC[2]>=5) || (BC[3]>=5))) || (REF=="G" && ((BC[0]>=5) || (BC[1]>=5) || (BC[3]>=5))) || (REF=="T" && ((BC[0]>=5) || (BC[1]>=5) || (BC[2]>=5))))' {output.vcf_file_gzip} > {output.vcf_SNPs}
        
        # indels won't be included, so save a separate file of indels, where IC >= 5 or DC >= 5
        bcftools filter -i "(strlen(REF) != strlen(ALT)) & (IC >= 5 | DC >= 5)" {output.vcf_file_gzip} > {output.vcf_indels}
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
        # convert the full VCF file to a BCF file to get only the lineage-defining positions according to the Coll 2014 scheme
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
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    params:
        lineage_SNP_info = os.path.join(primary_directory, references_dir, "phylogeny", "Coll2014_SNPs_all.csv"),
        F2_metric_script = os.path.join(primary_directory, scripts_dir, "calculate_F2_metric.py"),
        output_dir = output_dir,        
    output:
        F2_metric_output = f"{sample_out_dir}/lineage/F2_Coll2014.txt",
        minor_allele_fractions_output = temp(f"{sample_out_dir}/lineage/minor_allele_fractions.csv"),
        
        vcf_file_gunzip = temp(f"{sample_out_dir}/lineage/{{sample_ID}}_full.vcf"),
        fast_lineage_caller_output = f"{sample_out_dir}/lineage/fast_lineage_caller_output.txt",
    shell:
        """
        python3 -u {params.F2_metric_script} -i {params.output_dir}/{wildcards.sample_ID} -o {output.F2_metric_output} -O {output.minor_allele_fractions_output} --lineage-file {params.lineage_SNP_info}

        rm {input.bcf_file} {input.bcf_index_file} {input.vcf_lineage_positions}
        
        # fast-lineage-caller won't work on gzipped files, so need to unzip it first. 
        # It doesn't even error when you pass in a gzipped file. It just returns nothing, making debugging difficult
        gunzip -c {input.vcf_file_gzip} > {output.vcf_file_gunzip}

        fast-lineage-caller {output.vcf_file_gunzip} --pass --out {output.fast_lineage_caller_output}
        """