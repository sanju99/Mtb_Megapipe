import os, glob
import numpy as np
import pandas as pd

# define some paths to make the path names more readable
sample_out_dir = f"{output_dir}/{{sample_ID}}"
run_out_dir = f"{output_dir}/{{sample_ID}}/{{run_ID}}"

scripts_dir = config["scripts_dir"]
references_dir = config["references_dir"]

primary_directory = "/home/sak0914/Mtb_Megapipe"


rule convert_BAM_to_CRAM:
    input:
        bam_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.bam",
    output:
        cram_file = f"{sample_out_dir}/bam/{{sample_ID}}.dedup.cram",
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
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

        # delete the original BAM file
        rm {input.bam_file}
        """
        
        
        
rule convert_BAM_to_CRAM_subdirs:
    input:
        bam_file = f"{run_out_dir}/bam/{{run_ID}}.dedup.bam",
    output:
        cram_file = f"{run_out_dir}/bam/{{run_ID}}.dedup.cram",
    conda:
        f"{primary_directory}/envs/read_processing_aln.yaml"
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

        # delete the original BAM file
        rm {input.bam_file}
        """