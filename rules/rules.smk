import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"
run_out_dir = f"{output_dir}/{{sample_ID}}/{{run_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]

primary_directory = "/home/sak0914/Mtb_Megapipe"


rule get_input_FASTQ_files:
    group: 
        "sequential"
    output:
        fastq1 = f"{run_out_dir}/{{run_ID}}_1.fastq",
        fastq2 = f"{run_out_dir}/{{run_ID}}_2.fastq",
    params:
        sample_out_dir = sample_out_dir,
        fastq_dir = config["fastq_dir"],
        download_script = f"{primary_directory}/scripts/download_FASTQ.sh",
    run:        
        if download_public_FASTQ_dict[wildcards.sample_ID] == 1:
            shell("""
                module load sratoolkit/3.2.0

                # the script performs the same QC as in the next block
                bash {params.download_script} {params.sample_out_dir} {wildcards.run_ID}
            """)
        elif download_public_FASTQ_dict[wildcards.sample_ID] == 0:
            shell("""
                # copy the FASTQ files from the directory specified in the config file to the sample directory
                # they will be deleted in the next rule after performing adapter trimming, so they won't be doubly stored
                gunzip -c {params.fastq_dir}/{wildcards.run_ID}/{wildcards.run_ID}_R1.fastq.gz > {output.fastq1}
                gunzip -c {params.fastq_dir}/{wildcards.run_ID}/{wildcards.run_ID}_R2.fastq.gz > {output.fastq2}

                # first check that the original FASTQ files have the same numbers of lines
                FQ1_line_count=$(wc -l {output.fastq1} | awk '{{print $1}}')
                FQ2_line_count=$(wc -l {output.fastq2} | awk '{{print $1}}')
                
                # check that neither FASTQ file has no reads
                if [ $FQ1_line_count -eq 0 ] || [ $FQ2_line_count -eq 0 ]; then
                    echo "Error: At least one of the FASTQ files for $sample_ID/$run_ID has no reads"
                    exit 1
                # Compare the counts and raise an error if they are not equal 
                elif [ "$FQ1_line_count" -ne "$FQ2_line_count" ]; then
                    echo "Error: FASTQ files for $sample_ID/$run_ID have different line counts: $FQ1_line_count and $FQ2_line_count"
                    exit 1
                fi
                
                # compare paired end read files. If they are the same, then add to error list. Suppress output with -s tag, so it doesn't print out the differences
                # If the files are identical, the exit status is 0, and the condition is considered true, so an error will be returned.
                if cmp -s {output.fastq1} {output.fastq2}; then
                   echo "Error: {output.fastq1} and {output.fastq2} are duplicates"
                   exit 1
                fi
            """)


rule trim_adapters:
    input:
        fastq1 = f"{run_out_dir}/{{run_ID}}_1.fastq",
        fastq2 = f"{run_out_dir}/{{run_ID}}_2.fastq",
    output:
        fastq1_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R2.trimmed.fastq.gz",
        fastp_html = f"{run_out_dir}/fastp/fastp.html",
        fastp_json = f"{run_out_dir}/fastp/fastp.json"
    conda:
        f"{primary_directory}/envs/read_processing_aln_bwa.yaml"
    params:
        min_read_length = config["min_read_length"]
    shell:
        """
        fastp -i {input.fastq1} -I {input.fastq2} -o {output.fastq1_trimmed} -O {output.fastq2_trimmed} -h {output.fastp_html} -j {output.fastp_json} --length_required {params.min_read_length} --dedup --thread 8

        rm {input.fastq1} {input.fastq2}
        """



rule kraken_classification:
    input:
        fastq1_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R2.trimmed.fastq.gz",
    output:
        fastq1_trimmed_classified_unzipped = temp(f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq"),
        fastq2_trimmed_classified_unzipped = temp(f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq"),
        fastq1_trimmed_classified = f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq.gz",
        fastq2_trimmed_classified = f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq.gz",
        kraken_report = f"{run_out_dir}/kraken/kraken_report",
        kraken_classifications = temp(f"{run_out_dir}/kraken/kraken_classifications"),
    conda:
        f"{primary_directory}/envs/read_processing_aln_bwa.yaml"
    params:
        kraken_db = {config['kraken_db']},
        output_dir = output_dir,
        classified_out_string = f"{run_out_dir}/kraken/{{run_ID}}#.kraken.filtered.fastq"
    shell:
        """
        # --confidence is the minimum fraction of k-mers in a read that must match a given taxon for that read to be assigned to that taxon
        kraken2 --db {params.kraken_db} \
                --threads 8 \
                --confidence 0 \
                --paired {input.fastq1_trimmed} {input.fastq2_trimmed} \
                --gzip-compressed \
                --report {output.kraken_report} \
                --classified-out {params.classified_out_string} \
                --output {output.kraken_classifications}
                
        gzip -c {output.fastq1_trimmed_classified_unzipped} > {output.fastq1_trimmed_classified}
        gzip -c {output.fastq2_trimmed_classified_unzipped} > {output.fastq2_trimmed_classified}
        
        rm {input.fastq1_trimmed} {input.fastq2_trimmed}
        """
        

def compute_kraken_unclassified_percent(output_dir, sample_ID, run_ID):
    """
    This function determines whether or not a WGS run passes kraken classification. It returns True if the proportion of unclassified reads is less than or equal to the user-set maximum and False if the unclassified proportion exceeds it. Use percentage (out of 100).
    """

    kraken_report = pd.read_csv(f"{output_dir}/{sample_ID}/{run_ID}/kraken/kraken_report", sep='\t', header=None)

    # column 0 is the percentages, column 5 is the taxonomies. If unclassified is not in the column, there are no unclassified reads
    if 'unclassified' in kraken_report[5].values:
        kraken_unclassified_perc = kraken_report.loc[kraken_report[5]=='unclassified'][0].values[0]
    else:
        kraken_unclassified_perc = 0

    out_fName = f"{output_dir}/{sample_ID}/{run_ID}/kraken/unclassified_percent.txt"

    # if the file passes, create an output file
    if kraken_unclassified_perc <= config["kraken_unclassified_max"]:
        with open(out_fName, "w+") as file:
            file.write(str(kraken_unclassified_perc))
            file.write("\n")

    else:
        os.remove(f"{output_dir}/{sample_ID}/{run_ID}/kraken//{run_ID}_1.kraken.filtered.fastq.gz")
        os.remove(f"{output_dir}/{sample_ID}/{run_ID}/kraken//{run_ID}_2.kraken.filtered.fastq.gz")

    ## if not, do not create the file and delete the FASTQ files of classified reads because they are not needed
    #else:
    #    # write an empty file, appending so that it doesn't modify any existing content
    #    with open(out_fName, 'a'):
    #        # Update the timestamp to simulate the behavior of touch in bash
    #        os.utime(out_fName, None)


checkpoint check_kraken_unclassified:
    input:
        kraken_report = f"{run_out_dir}/kraken/kraken_report",
    output: 
        kraken_pass_file = f"{run_out_dir}/kraken/unclassified_percent.txt"
    params:
        output_dir = output_dir,
    run:
        compute_kraken_unclassified_percent(params.output_dir, wildcards.sample_ID, wildcards.run_ID)

        # Check which samples passed the threshold and record them
        if os.path.isfile(output.kraken_pass_file):
            # If the file exists and has content, it's a valid sample
            print(f"Sample has less than or equal to {config['kraken_unclassified_max']}% unclassified reads")
        else:
            print(f"Sample has more than {config['kraken_unclassified_max']}% unclassified reads. Halting execution")



rule fastlin_typing:
    input:
        kraken_pass_file = f"{run_out_dir}/kraken/unclassified_percent.txt",
        fastq1_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq.gz",
        fastq2_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq.gz",
    output:
        fastq1_trimmed_classified_gzipped = temp(f"{run_out_dir}/fastlin/{{run_ID}}_1.fastq.gz"),
        fastq2_trimmed_classified_gzipped = temp(f"{run_out_dir}/fastlin/{{run_ID}}_2.fastq.gz"),
        fastlin_dir = directory(f"{run_out_dir}/fastlin"),
        fastlin_output = f"{run_out_dir}/fastlin/output.txt"
    conda:
        f"{primary_directory}/envs/read_processing_aln_bwa.yaml"
    params:
        fastlin_barcodes = os.path.join(primary_directory, references_dir, "phylogeny", "MTBC_barcodes.tsv"),
    shell:
        """
        cp {input.fastq1_trimmed_classified} {output.fastq1_trimmed_classified_gzipped}
        cp {input.fastq2_trimmed_classified} {output.fastq2_trimmed_classified_gzipped}
        
        # there's one TRUST sample with coverage > 900
        fastlin -d {output.fastlin_dir} -b {params.fastlin_barcodes} -o {output.fastlin_output} -x 100
        """


def does_sample_pass_fastlin(output_dir, sample_ID):
    """
    Use this function only for cases when there are multiple distinct sequencing runs for a single sample. Before merging the individual BAMs, check that they are assigned the same lineage according to fastlin. This is to mitigate the risk of potentially mislabeled WGS runs being merged together. 
    """    
    fastlin_outputs = pd.concat([pd.read_csv(fName, sep='\t') for fName in glob.glob(f"{output_dir}/{sample_ID}/*/fastlin/output.txt")]).reset_index(drop=True)

    if len(fastlin_outputs) == 0:
        raise ValueError(f"There are no fastlin outputs for {output_dir}/{sample_ID}")

    # check the coverage in the k_cov column. This can save time by identifying no-coverage samples and halting the pipeline before bwa-mem
    # these are probably targeted sequencing samples, so the genome-wide coverage is estimated at 0 by fastlin
    for i, row in fastlin_outputs.iterrows():

        if row['k_cov'] == 0:
            low_coverage_runs += 1
            print(f"{sample_ID}/{row['#sample']} has a k-mer coverage of {row['k_cov']}")

            # delete the kraken-filtered FASTQ files to prevent bwa-mem from being run
            os.remove(f"{output_dir}/{sample_ID}/{row['#sample']}/kraken/{row['#sample']}_1.kraken.filtered.fastq.gz")
            os.remove(f"{output_dir}/{sample_ID}/{row['#sample']}/kraken/{row['#sample']}_2.kraken.filtered.fastq.gz")

        # split lineage from median k-mer occurrence
        elif ',' not in row['lineages']:

            # most common lineage is the same as the only lineage present
            fastlin_outputs.loc[i, ['lineage', 'lineage_koccur', 'most_common_lineage', 'most_common_lineage_koccur']] = [row['lineages'].split(' ')[0], row['lineages'].split(' ')[1].replace('(', '').replace(')', ''), row['lineages'].split(' ')[0], row['lineages'].split(' ')[1].replace('(', '').replace(')', '')]
        
        else:
    
            fastlin_lineage_lst = []
            fastlin_median_occur_lst = []
            
            for single_lineage in row['lineages'].split(', '):
                fastlin_lineage_lst.append(single_lineage.split(' ')[0])
                fastlin_median_occur_lst.append(single_lineage.split(' ')[1].replace('(', '').replace(')', ''))

            # sort by alpha-numeric order
            # need indices to sort the k-mer occurrences too
            alpha_sorted_idx = np.argsort(fastlin_lineage_lst)
            alpha_sorted_koccur_lst = [fastlin_median_occur_lst[idx] for idx in alpha_sorted_idx]

            # convert values in fastlin_median_occur_lst from strings to ints
            fastlin_median_occur_lst = np.array(fastlin_median_occur_lst).astype(int)

            # also add a column for the more common lineage in cases where there are multiple lineages assigned
            most_common_lineage = fastlin_lineage_lst[np.argmax(fastlin_median_occur_lst)]
            
            fastlin_outputs.loc[i, ['lineage', 'lineage_koccur', 'most_common_lineage', 'most_common_lineage_koccur']] = [','.join(np.sort(fastlin_lineage_lst)), ','.join(alpha_sorted_koccur_lst), most_common_lineage, np.max(fastlin_median_occur_lst)]

    fastlin_outputs['most_common_lineage_koccur'] = fastlin_outputs['most_common_lineage_koccur'].astype(int)
    
    # if the most common lineages across the runs match, then it probably indicates low-level contamination by another lineage, not a lineage mixture. We want to keep these
    if fastlin_outputs['most_common_lineage'].nunique() > 1:
        # return False
        print(f"Halting pipeline for {wildcards.sample_ID} because the different WGS runs for it have different lineages assigned by fastlin")
        exit()



checkpoint check_fastlin:
    input:
        fastlin_output = f"{run_out_dir}/fastlin/output.txt"
    params:
        output_dir = output_dir,
    run:
        does_sample_pass_fastlin(params.output_dir, wildcards.sample_ID)



def get_primary_lineage_from_fastlin_output(fastlin_fName):

    try:
        df = pd.read_csv(fastlin_fName, sep='\t')
    except:
        return None

    # check that everything is not NA. This can happen if there are too few reads in the FASTQ files to get anything
    if len(df.dropna(axis=1, how='all')) == 0:
        return None

    # check that there is a single sample and it is paired reads
    assert len(df) == 1
    assert df['data_type'].values[0] == 'paired' 

    for i, row in df.iterrows():
        if ',' not in row['lineages']:
            df.loc[i, ['lineage', 'lineage_koccur']] = [row['lineages'].split(' ')[0], row['lineages'].split(' ')[1].replace('(', '').replace(')', '')]
        else:
    
            fastlin_lineage_lst = []
            fastlin_median_occur_lst = []
            
            for single_lineage in row['lineages'].split(', '):
                fastlin_lineage_lst.append(single_lineage.split(' ')[0])
                fastlin_median_occur_lst.append(single_lineage.split(' ')[1].replace('(', '').replace(')', ''))
    
            # need indices to sort the k-mer occurrences too
            sorted_idx = np.argsort(fastlin_lineage_lst)
            sorted_koccur_lst = [fastlin_median_occur_lst[idx] for idx in sorted_idx]
            
            df.loc[i, ['lineage', 'lineage_koccur']] = [','.join(np.sort(fastlin_lineage_lst)), ','.join(sorted_koccur_lst)] 

    full_lineage = df['lineage'].values[0]

    if full_lineage.isnumeric():
        # take the first value (i.e. take 2 if the lineage is 2.2.1)
        primary_lineage = full_lineage[0]
    else:
        primary_lineage = full_lineage

    return primary_lineage



rule align_reads_mark_duplicates:
    input:
        # require the fastlin output file as an input so that the fastlin rule gets run
        fastlin_output = f"{run_out_dir}/fastlin/output.txt",
        fastq1_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq.gz",
        fastq2_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq.gz",
    output:
        sam_file = temp(f"{run_out_dir}/bam/{{run_ID}}.sam"),
        bam_file = temp(f"{run_out_dir}/bam/{{run_ID}}.bam"),
        bam_index_file = temp(f"{run_out_dir}/bam/{{run_ID}}.bam.bai"),
        bam_file_markDup = f"{run_out_dir}/bam/{{run_ID}}.markDup.bam",
        bam_file_markDup_metrics = f"{run_out_dir}/bam/{{run_ID}}.markDup.bam.metrics",
        bam_index_file_markDup = f"{run_out_dir}/bam/{{run_ID}}.markDup.bam.bai",
    params:
        output_dir = output_dir,
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
    conda:
        # f"{primary_directory}/envs/read_processing_aln_bwa.yaml"
        f"{primary_directory}/envs/read_processing_aln.yaml"
    shell:
        """
        # align reads to the reference genome sequence. The RG name specifies the read group name, which is necessary if you are merging multiple WGS runs into a single BAM file
        bwa-mem2 mem -M -R "@RG\\tID:{wildcards.run_ID}\\tSM:{wildcards.sample_ID}" -t 8 {params.ref_genome} {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} > {output.sam_file}
        
        # bwa mem -M -R "@RG\\tID:{wildcards.run_ID}\\tSM:{wildcards.sample_ID}" -t 8 {params.ref_genome} {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} > {output.sam_file}

        # sort alignment and convert to bam file
        samtools view -b {output.sam_file} | samtools sort > {output.bam_file}

        # index alignment, which creates a .bai index file
        samtools index {output.bam_file}

        # -Xmx6g specifies to allocate 6 GB. Don't remove duplicate reads, just mark them
        picard -Xmx10g MarkDuplicates I={output.bam_file} O={output.bam_file_markDup} REMOVE_DUPLICATES=false M={output.bam_file_markDup_metrics} ASSUME_SORT_ORDER=coordinate READ_NAME_REGEX='(?:.*.)?([0-9]+)[^.]*.([0-9]+)[^.]*.([0-9]+)[^.]*$'

        # index the alignment with samtools
        samtools index {output.bam_file_markDup}
        
        rm {input.fastq1_trimmed_classified}
        rm {input.fastq2_trimmed_classified}
        """


rule get_BAM_file_depths:
    input:
        # single_run_BAM_file = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.markDup.bam" for run_ID in sample_run_dict[wildcards.sample_ID] if os.path.isfile(f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.markDup.bam")],
        single_run_BAM_file = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.markDup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]]
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_out_dir = sample_out_dir,
    output:
        depth_file = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv.gz",
    conda:
        f"{primary_directory}/envs/read_processing_aln_bwa.yaml"
    shell:
        """
        # get all runs associated with this sample_ID and compute depth
        # -a computes depth at all positions, not just those with non-zero depth
        # -Q is for minimum mapping quality: use 1, so that multiply mapped reads aren't counted. These have mapping quality of 0
        # By default, reads that have any of the flags UNMAP, SECONDARY, QCFAIL, or DUP set are skipped in samtools depth v1.21
        samtools depth -a -Q 1 {input.single_run_BAM_file} | gzip -c > {output.depth_file}
        """


rule get_BAMs_passing_QC_thresholds:
    input:
        depth_file_gzip = f"{sample_out_dir}/bam/{{sample_ID}}.depth.tsv.gz", # contains depths for all BAM files for all WGS runs
        # single_run_BAM_file = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.markDup.bam" for run_ID in sample_run_dict[wildcards.sample_ID] if os.path.isfile(f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.markDup.bam")],
        single_run_BAM_file = lambda wildcards: [f"{output_dir}/{wildcards.sample_ID}/{run_ID}/bam/{run_ID}.markDup.bam" for run_ID in sample_run_dict[wildcards.sample_ID]]
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
        python3 -u {params.BAM_depth_QC_script} -i {input.depth_file_gzip} -b {input.single_run_BAM_file} --median-depth {params.median_depth} --min-cov {params.min_cov} --genome-cov-prop {params.genome_cov_prop}
        """



rule merge_BAMs:
    input:
        pass_BAMs_file = f"{sample_out_dir}/bam/pass_BAMs.txt",
    output:
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam",
        merged_bam_index_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam.bai",
    conda:
        f"{primary_directory}/envs/read_processing_aln_bwa.yaml"
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
        merged_bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam",
    output:
        vcf_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.vcf"),
        vcf_file_full = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
        fasta_file = temp(f"{sample_out_dir}/pilon/{{sample_ID}}.fasta"),
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    params:
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
        sample_pilon_dir = f"{sample_out_dir}/pilon",
    conda:
         f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        pilon -Xmx10g --minmq 1 --genome {params.ref_genome} --bam {input.merged_bam_file} --output {wildcards.sample_ID} --outdir {params.sample_pilon_dir} --variant
            
        # left-align indels and drop duplicates, then gzip the full VCF file 
        # this affects those cases where the position of the indel is ambiguous
        # however, because of the shifting positions, the position of the indel can change, so need to sort it
        bcftools norm --rm-dup none --fasta-ref {params.ref_genome} {output.vcf_file} | bcftools sort | gzip -c > {output.vcf_file_full}
        
        bcftools view --types snps,mnps,indels,other {output.vcf_file_full} > {output.vcf_file_variants_only}
        """



rule pilon_get_variants_only_file:
    input:
        vcf_file_full = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    output:
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    conda:
        f"{primary_directory}/envs/variant_calling.yaml"
    threads:
        1
    shell:
        """        
        bcftools view --types snps,mnps,indels,other {input.vcf_file_full} > {output.vcf_file_variants_only}
        """


rule lineage_typing:
    input:
        vcf_file_full = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf",     
    output:
        fast_lineage_caller_output = f"{sample_out_dir}/lineage/fast_lineage_caller_output.txt",
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    params:
        lineage_SNP_info = os.path.join(primary_directory, references_dir, "phylogeny", "Coll2014_SNPs_all.csv"),
        output_dir = output_dir,   
    shell:
        """        
        fast-lineage-caller {input.vcf_file_full} --pass --out {output.fast_lineage_caller_output}
        
        # then gzip the full vcf file
        gzip {input.vcf_file_full}
        """
        
        

rule create_F2_helper_files:
    input:
        fast_lineage_caller_output = f"{sample_out_dir}/lineage/fast_lineage_caller_output.txt",
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
    params:
        output_dir = output_dir,
    output:
        bcf_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf",
        bcf_index_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf.csi",
        vcf_lineage_positions = f"{sample_out_dir}/lineage/{{sample_ID}}_lineage_positions.vcf",
    conda:
         f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        chrom_name=$(zgrep -v '^#' {input.vcf_file_gzip} | awk 'NR==1 {{print $1}}')
    
        if [ "$chrom_name" = "Chromosome" ]; then
            lineage_pos_file="/home/sak0914/Mtb_Megapipe/references/phylogeny/Coll2014_positions_all.txt"
        else
            lineage_pos_file="/home/sak0914/Mtb_Megapipe/references/phylogeny/Coll2014_positions_all_refseq.txt"
        fi
        
        # convert the full VCF file to a BCF fileto get only the lineage-defining positions according to the Coll 2014 scheme
        bcftools view {input.vcf_file_gzip} -O b -o {output.bcf_file}

        # index bcf file
        bcftools index {output.bcf_file}

        # create VCF file of just the lineage positions, which will be used by the F2 metric script. Per the documentation, if --regions-file is a tab-delimited file, then it needs two columns (CHROM and POS), and POS is 1-indexed and inclusive
        # THIS IS DIFFERENT BEHAVIOR FROM IF IT WAS A BED FILE OR IF YOU USE BEDTOOLS. IN BOTH OF THOSE CASES, YOU NEED THREE COLUMNS (CHROM, BEG, AND END), AND THEY ARE 0-INDEXED WITH END BEING EXCLUSIVE (I.E. HALF-OPEN)
        bcftools view {output.bcf_file} --regions-file $lineage_pos_file -O v -o {output.vcf_lineage_positions}   
        """



rule compute_F2_metric:
    input:
        bcf_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf",
        bcf_index_file = f"{sample_out_dir}/lineage/{{sample_ID}}.bcf.csi",
        vcf_lineage_positions = f"{sample_out_dir}/lineage/{{sample_ID}}_lineage_positions.vcf",
    output:
        F2_metric_output = f"{sample_out_dir}/lineage/F2_Coll2014.txt",
        minor_allele_fractions_output = temp(f"{sample_out_dir}/lineage/minor_allele_fractions.csv"),
    params:
        lineage_SNP_info = os.path.join(primary_directory, references_dir, "phylogeny", "Coll2014_SNPs_all.csv"),
        # F2_metric_script = os.path.join(primary_directory, scripts_dir, "calculate_strain_mixing_VCF.py"),
        F2_metric_script = os.path.join(primary_directory, scripts_dir, "calculate_F2_metric.py"),
        output_dir = output_dir,
    shell:
        """
        python3 -u {params.F2_metric_script} -i {params.output_dir}/{wildcards.sample_ID} -o {output.F2_metric_output} -O {output.minor_allele_fractions_output} --lineage-file {params.lineage_SNP_info} #--num_mixed_lineages 2

        rm {input.bcf_file} {input.bcf_index_file} {input.vcf_lineage_positions}
        """



rule exclude_low_confidence_sites:
    input:
        vcf_file_variants_only = f"{sample_out_dir}/pilon/{{sample_ID}}_variants.vcf",
    output:
        vcf_file_exclude_low_quality = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_excludeLowConf.vcf",
    params:
        exclude_regions = os.path.join(primary_directory, references_dir, "BED_files/exclude_regions_50bp.bed"),
    threads:
        1
    conda:
         f"{primary_directory}/envs/variant_calling.yaml"
    shell:
        """
        # source activate /home/sak0914/anaconda3/envs/python_bioinformatics_utils
        
        chrom_name=$(zgrep -v '^#' {input.vcf_file_variants_only} | awk 'NR==1 {{print $1}}')
    
        if [ "$chrom_name" = "Chromosome" ]; then
            exclude_regions="/home/sak0914/Mtb_Megapipe/references/BED_files/exclude_regions_50bp.bed"
        else
            exclude_regions="/home/sak0914/Mtb_Megapipe/references/BED_files/exclude_regions_50bp_refseq.bed"
        fi
        
        bcftools view -T ^$exclude_regions {input.vcf_file_variants_only} > {output.vcf_file_exclude_low_quality}
        """
        
        
        

rule combine_codon_variants:
    input:
        vcf_file_exclude_low_quality = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_excludeLowConf.vcf",
    output:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_excludeLowConf_combinedCodons.vcf",
    params:
        combine_codon_variants_script = os.path.join(primary_directory, scripts_dir, "combine_codon_variants.py"),
    threads:
        1
    shell:
        """
        # this creates _variants_combined_codons.vcf in the same directory as the input VCF. You can also specify the -o flag if you want to write the output file somewhere else
        python3 -u {params.combine_codon_variants_script} -i {input.vcf_file_exclude_low_quality}        
        """



rule annotate_variants_snpEff:
    input:
        vcf_file_variants_combinedCodons = f"{sample_out_dir}/pilon/{{sample_ID}}_variants_excludeLowConf_combinedCodons.vcf",
    output:
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_excludeLowConf_combinedCodons.eff.vcf",
    conda:
        f"{primary_directory}/envs/variant_annotation.yaml"
    params:
        snpEff_db = config['snpEff_db'],
    threads:
        1
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
        vcf_file_variants_combinedCodons_annot = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants_excludeLowConf_combinedCodons.eff.vcf",
    output:
        variants_file_tsv = f"{sample_out_dir}/WHO_resistance/{{sample_ID}}_variants.tsv",
    conda:
        f"{primary_directory}/envs/variant_annotation.yaml"
    threads:
        1
    shell:
        """
        chrom_name=$(zgrep -v '^#' {input.vcf_file_variants_combinedCodons_annot} | awk 'NR==1 {{print $1}}')
    
        if [ "$chrom_name" = "Chromosome" ]; then
            WHO_resistance_bed_file="/home/sak0914/Mtb_Megapipe/references/WHO_catalog_resistance/regions.bed"
            exclude_regions_bed_file="/home/sak0914/Mtb_Megapipe/references/BED_files/RLC_Regions.Plus.LowPmapK50E4.H37Rv.bed"
        else
            WHO_resistance_bed_file="/home/sak0914/Mtb_Megapipe/references/WHO_catalog_resistance/regions_refseq.bed"
            exclude_regions_bed_file="/home/sak0914/Mtb_Megapipe/references/BED_files/RLC_Regions.Plus.LowPmapK50E4.H37Rv.RefSeq.bed"
        fi

        # exclude regions is regions with low EBR and low 50-mer mappability
        # bedtools subtract -a {input.vcf_file_variants_combinedCodons_annot} -b $exclude_regions_bed_file -header | bedtools intersect -a '-' -b $WHO_resistance_bed_file -header | SnpSift extractFields '-' POS REF ALT FILTER QUAL IMPRECISE AF DP BQ MQ IC DC ANN -e "" > {output.variants_file_tsv}
        
        bedtools intersect -a {input.vcf_file_variants_combinedCodons_annot} -b $WHO_resistance_bed_file -header | SnpSift extractFields '-' POS REF ALT FILTER QUAL IMPRECISE AF DP BQ MQ IC DC ANN -e "" > {output.variants_file_tsv}
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
    threads:
        1
    shell:
        """
        python3 -u {params.process_variants_WHO_catalog_script} -i {input.variants_file_tsv} -o {output.variants_file_tsv_annot}
        rm {input.variants_file_tsv}

        # get resistance predictions -- any Group 1 or 2 variant that passes QC leads to a prediction of R for a given drug. If not, predicted S
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o {params.output_file_basename}
        python3 -u {params.get_WHO_resistance_predictions_script} -i {output.variants_file_tsv_annot} -o {params.output_file_basename} --AF-thresh 0.25
        """
        
        
        
rule convert_BAM_to_CRAM:
    input:
        vcf_file_gzip = f"{sample_out_dir}/pilon/{{sample_ID}}_full.vcf.gz",
        bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.bam",
    output:
        cram_file = f"{sample_out_dir}/bam/{{sample_ID}}.cram",
    conda:
        f"{primary_directory}/envs/read_processing_aln_bwa.yaml"
    shell:
        """        
        # determine the appropriate FASTA file by getting the name of the reference chromosome used for alignment
        chrom_name=$(samtools idxstats {input.bam_file}  | cut -f1 | head -1)

        if [ "$chrom_name" = "Chromosome" ]; then
            ref_fasta="/home/sak0914/Mtb_Megapipe/references/ref_genome/H37Rv_NC_000962.3.fna"
        elif [ "$chrom_name" = "NC_000962.3" ]; then
            ref_fasta="/home/sak0914/Mtb_Megapipe/references/ref_genome/refseq.fna"
        else
            echo "CHROM name $chrom_name is invalid"
        fi
        
        # create CRAM file
        samtools view -T $ref_fasta -C -o {output.cram_file} {input.bam_file}

        # delete the original BAM file and index file. The index file has to be regenerated when the new BAM file is created as well
        rm {input.bam_file} {input.bam_file}.bai
        """