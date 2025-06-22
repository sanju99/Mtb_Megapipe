#!/bin/bash
#SBATCH -c 1
#SBATCH -t 2-06:00
#SBATCH -p medium
#SBATCH --mem=5G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

source activate bioinformatics

set -o errexit # any error will cause the shell script to exit immediately. This is not native bash behavior

if ! [ $# -eq 1 ]; then
    echo "Please pass in 1 command line argument: a text file of the full names of BAM files to compress to CRAM"
    exit
fi

input_file=$1

# Read the TSV file line by line, skiping the header. IFS sets the field separator. Here, it is tab
while IFS=$'\t', read -r BAM_fName
do

    if [ -f $BAM_fName ]; then

        # determine the appropriate FASTA file by getting the name of the reference chromosome used for alignment
        chrom_name=$(samtools idxstats $BAM_fName | cut -f1 | head -1)

        if [ "$chrom_name" = "Chromosome" ]; then
            ref_fasta="/home/sak0914/Mtb_Megapipe/references/ref_genome/H37Rv_NC_000962.3.fna"
        elif [ "$chrom_name" = "NC_000962.3" ]; then
            ref_fasta="/home/sak0914/Mtb_Megapipe/references/ref_genome/refseq.fna"
        else
            echo "CHROM name $chrom_name is invalid"
        fi

        CRAM_fName="${BAM_fName%.bam}.cram"

        # create CRAM file
        samtools view -T $ref_fasta -C -o $CRAM_fName $BAM_fName

        # delete the original BAM file
        rm $BAM_fName

        echo "Finished converting $BAM_fName"
        
    else
        echo "Already finished $BAM_fName"
        
    fi

done < "$input_file"