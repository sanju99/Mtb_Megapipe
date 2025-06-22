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
        fastq1 = f"{run_out_dir}/{{run_ID}}_R1.fastq.gz",
        fastq2 = f"{run_out_dir}/{{run_ID}}_R2.fastq.gz",

        fastq1_unzipped = temp(f"{run_out_dir}/{{run_ID}}_1.fastq"),
        fastq2_unzipped = temp(f"{run_out_dir}/{{run_ID}}_2.fastq"),
    params:
        sample_out_dir = sample_out_dir,
        fastq_dir = config["fastq_dir"],
        download_script = f"{primary_directory}/scripts/download_FASTQ.sh",
    run:        
        if download_public_FASTQ_dict[wildcards.sample_ID] == 1:
            shell("""
                module load sratoolkit/2.10.7

                # the script performs the same QC as in the next block
                bash {params.download_script} {params.sample_out_dir} {wildcards.run_ID}
            """)
        elif download_public_FASTQ_dict[wildcards.sample_ID] == 0:
            shell("""
                # copy the FASTQ files from the directory specified in the config file to the sample directory
                # they will be deleted in the next rule after performing adapter trimming, so they won't be doubly stored
                cp {params.fastq_dir}/{wildcards.run_ID}/{wildcards.run_ID}_R1.fastq.gz {output.fastq1}
                cp {params.fastq_dir}/{wildcards.run_ID}/{wildcards.run_ID}_R2.fastq.gz {output.fastq2}

                gunzip -c {output.fastq1} > {output.fastq1_unzipped}
                gunzip -c {output.fastq2} > {output.fastq2_unzipped}

                # first check that the original FASTQ files have the same numbers of lines
                FQ1_line_count=$(wc -l {output.fastq1_unzipped} | awk '{{print $1}}')
                FQ2_line_count=$(wc -l {output.fastq2_unzipped} | awk '{{print $1}}')
                
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
                if cmp -s {output.fastq1_unzipped} {output.fastq2_unzipped}; then
                   echo "Error: {output.fastq1_unzipped} and {output.fastq2_unzipped} are duplicates"
                   exit 1
                fi
            """)


rule trim_adapters:
    input:
        fastq1 = f"{run_out_dir}/{{run_ID}}_R1.fastq.gz",
        fastq2 = f"{run_out_dir}/{{run_ID}}_R2.fastq.gz",
    output:
        fastq1_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R2.trimmed.fastq.gz",
        fastp_html = f"{run_out_dir}/fastp/fastp.html",
        fastp_json = f"{run_out_dir}/fastp/fastp.json"
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    params:
        min_read_length = config["min_read_length"]
    shell:
        """
        fastp -i {input.fastq1} -I {input.fastq2} -o {output.fastq1_trimmed} -O {output.fastq2_trimmed} -h {output.fastp_html} -j {output.fastp_json} --length_required {params.min_read_length} --dedup --thread 8

        # rm {input.fastq1} {input.fastq2}
        """



rule kraken_classification:
    input:
        fastq1_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R1.trimmed.fastq.gz",
        fastq2_trimmed = f"{run_out_dir}/fastp/{{run_ID}}.R2.trimmed.fastq.gz",
    output:
        fastq1_trimmed_classified = f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq",
        fastq2_trimmed_classified = f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq",
        kraken_report = f"{run_out_dir}/kraken/kraken_report",
        kraken_classifications = temp(f"{run_out_dir}/kraken/kraken_classifications"),
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
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
        
        # rm {input.fastq1_trimmed} {input.fastq2_trimmed}
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
        os.remove(f"{output_dir}/{sample_ID}/{run_ID}/kraken//{run_ID}_1.kraken.filtered.fastq")
        os.remove(f"{output_dir}/{sample_ID}/{run_ID}/kraken//{run_ID}_2.kraken.filtered.fastq")

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
        # this is how you access the checkpoint outputs. Each checkpoint is stored in the checkpoints global variable
        kraken_pass_file = f"{run_out_dir}/kraken/unclassified_percent.txt",
        fastq1_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq",
        fastq2_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq",
    output:
        fastq1_trimmed_classified_gzipped = temp(f"{run_out_dir}/fastlin/{{run_ID}}_1.fastq.gz"),
        fastq2_trimmed_classified_gzipped = temp(f"{run_out_dir}/fastlin/{{run_ID}}_2.fastq.gz"),
        fastlin_dir = directory(f"{run_out_dir}/fastlin"),
        fastlin_output = f"{run_out_dir}/fastlin/output.txt"
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    params:
        fastlin_barcodes = os.path.join(primary_directory, references_dir, "phylogeny", "MTBC_barcodes.tsv"),
    shell:
        """
        gzip -c {input.fastq1_trimmed_classified} > {output.fastq1_trimmed_classified_gzipped}
        gzip -c {input.fastq2_trimmed_classified} > {output.fastq2_trimmed_classified_gzipped}
        
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
            os.remove(f"{output_dir}/{sample_ID}/{row['#sample']}/kraken/{row['#sample']}_1.kraken.filtered.fastq")
            os.remove(f"{output_dir}/{sample_ID}/{row['#sample']}/kraken/{row['#sample']}_2.kraken.filtered.fastq")

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



rule align_reads_mark_duplicates_bwa_mem2:
    input:
        # require the fastlin output file as an input so that the fastlin rule gets run
        fastlin_output = f"{run_out_dir}/fastlin/output.txt",
        fastq1_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq",
        fastq2_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq",
    output:
        sam_file = temp(f"{run_out_dir}/bwa_mem2/{{run_ID}}.sam"),
        bam_file = temp(f"{run_out_dir}/bwa_mem2/{{run_ID}}.bam"),
        bam_index_file = temp(f"{run_out_dir}/bwa_mem2/{{run_ID}}.bam.bai"),
        bam_file_markDup = f"{run_out_dir}/bwa_mem2/{{run_ID}}.markDup.bam",
        bam_file_markDup_metrics = f"{run_out_dir}/bwa_mem2/{{run_ID}}.markDup.bam.metrics",
        bam_index_file_markDup = f"{run_out_dir}/bwa_mem2/{{run_ID}}.markDup.bam.bai",
    params:
        output_dir = output_dir,
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
    shell:
        """
        # index reference genome (which is required before aligning reads)
        # bwa-mem2 index {params.ref_genome}

        # align reads to the reference genome sequence. The RG name specifies the read group name, which is necessary if you are merging multiple WGS runs into a single BAM file
        bwa-mem2 mem -k 80 -M -R "@RG\\tID:{wildcards.run_ID}\\tSM:{wildcards.sample_ID}" -t 8 {params.ref_genome} {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} > {output.sam_file}

        # sort alignment and convert to bam file
        samtools view -b {output.sam_file} | samtools sort > {output.bam_file}

        # index alignment, which creates a .bai index file
        samtools index {output.bam_file}

        # -Xmx6g specifies to allocate 6 GB. Don't remove duplicate reads, just mark them
        picard -Xmx10g MarkDuplicates I={output.bam_file} O={output.bam_file_markDup} REMOVE_DUPLICATES=false M={output.bam_file_markDup_metrics} ASSUME_SORT_ORDER=coordinate READ_NAME_REGEX='(?:.*.)?([0-9]+)[^.]*.([0-9]+)[^.]*.([0-9]+)[^.]*$'

        # index the deduplicated alignment with samtools, which will create a dedup_bam_file.bai file
        samtools index {output.bam_file_markDup}
        """


rule align_reads_mark_duplicates_bwa_mem:
    input:
        # require the fastlin output file as an input so that the fastlin rule gets run
        fastlin_output = f"{run_out_dir}/fastlin/output.txt",
        fastq1_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_1.kraken.filtered.fastq",
        fastq2_trimmed_classified=f"{run_out_dir}/kraken/{{run_ID}}_2.kraken.filtered.fastq",
    output:
        sam_file = temp(f"{run_out_dir}/bwa_mem/{{run_ID}}.sam"),
        bam_file = temp(f"{run_out_dir}/bwa_mem/{{run_ID}}.bam"),
        bam_index_file = temp(f"{run_out_dir}/bwa_mem/{{run_ID}}.bam.bai"),
        bam_file_markDup = f"{run_out_dir}/bwa_mem/{{run_ID}}.markDup.bam",
        bam_file_markDup_metrics = f"{run_out_dir}/bwa_mem/{{run_ID}}.markDup.bam.metrics",
        bam_index_file_markDup = f"{run_out_dir}/bwa_mem/{{run_ID}}.markDup.bam.bai",
    params:
        output_dir = output_dir,
        ref_genome = os.path.join(primary_directory, references_dir, "ref_genome", "H37Rv_NC_000962.3.fna"),
    conda:
        f"{primary_directory}/envs/bwa_mem.yaml"
    shell:
        """
        # align reads to the reference genome sequence. The RG name specifies the read group name, which is necessary if you are merging multiple WGS runs into a single BAM file
        bwa mem -k 80 -M -R "@RG\\tID:{wildcards.run_ID}\\tSM:{wildcards.sample_ID}" -t 8 {params.ref_genome} {input.fastq1_trimmed_classified} {input.fastq2_trimmed_classified} > {output.sam_file}

        # sort alignment and convert to bam file
        samtools view -b {output.sam_file} | samtools sort > {output.bam_file}

        # index alignment, which creates a .bai index file
        samtools index {output.bam_file}

        # -Xmx6g specifies to allocate 6 GB. Don't remove duplicate reads, just mark them
        picard -Xmx10g MarkDuplicates I={output.bam_file} O={output.bam_file_markDup} REMOVE_DUPLICATES=false M={output.bam_file_markDup_metrics} ASSUME_SORT_ORDER=coordinate READ_NAME_REGEX='(?:.*.)?([0-9]+)[^.]*.([0-9]+)[^.]*.([0-9]+)[^.]*$'

        # index the deduplicated alignment with samtools, which will create a dedup_bam_file.bai file
        samtools index {output.bam_file_markDup}
        """