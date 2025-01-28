#!/bin/bash 
#SBATCH -c 1
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=2G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

set -o errexit # any error will cause the shell script to exit immediately. This is not native bash behavior
# source activate bioinformatics

# input_file="/home/sak0914/MtbQuantCNN/VCFs_remove_duplicate_indels.csv"
input_file=$1
out_dir="/n/data1/hms/dbmi/farhat/rollingDB/genomic_data"

# Check if the output directory exists. If not, raise an error
if [ ! -d "$out_dir" ]; then
    echo "Output directory $out_dir doesn't exist!"
    exit 1
fi

# Read the TSV file line by line, skiping the header. IFS sets the field separator. Here, it is tab
while IFS=$',', read -r sample_ID
do

    if [ ! -d "$out_dir/$sample_ID/WHO_resistance" ]; then
        mkdir "$out_dir/$sample_ID/WHO_resistance"
    fi

    if [ ! -f "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.vcf" ]; then
    
        python3 -u "/home/sak0914/Mtb_Megapipe/scripts/combine_codon_variants.py" -i "$out_dir/$sample_ID/pilon/${sample_ID}_variants.vcf" -o "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.vcf"

        echo "Finished $sample_ID"

    else

        echo "Already did $sample_ID"

    fi

done < "$input_file"