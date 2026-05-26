#!/bin/bash
#SBATCH -c 2
#SBATCH -t 2-00:00
#SBATCH -p medium
#SBATCH --mem=10G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

source activate /home/sak0914/Mtb_Megapipe/.snakemake/conda/08dcc4c031abc315bb6cdf6e29a2302a_

set -o errexit # any error will cause the shell script to exit immediately. This is not native bash behavior

if ! [ $# -eq 1 ]; then
    echo "Please pass in 1 command line argument: a text file of the full names of SAM files to compress to BAM"
    exit
fi

input_file=$1

# Read the TSV file line by line, skiping the header. IFS sets the field separator. Here, it is tab
while IFS=$'\t' read -r SAM_fName
do

    BAM_fName="${SAM_fName%.sam}.bam"

    if [ -f $SAM_fName ]; then

        # convert to BAM and sort
        samtools view -b $SAM_fName --threads 4 | samtools sort --threads 4 > $BAM_fName

        # index alignment, which creates a .bai index file
        samtools index $BAM_fName

        # delete the original BAM file and the index file. The index has to be remade from the new BAM file
        rm $SAM_fName

        echo "Finished converting $SAM_fName"
        
    else
        echo "Already finished $SAM_fName"
        
        # make BAM index if not there
        if [ ! -f "$BAM_fName.bai" ]; then
            
            samtools index $BAM_fName
            
        fi
        
    fi

done < "$input_file"