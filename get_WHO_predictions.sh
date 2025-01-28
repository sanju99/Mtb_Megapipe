#!/bin/bash 
#SBATCH -c 1
#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mem=5G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

set -o errexit # any error will cause the shell script to exit immediately. This is not native bash behavior
source activate bioinformatics

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

    if [ ! -f "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_pred_AF_thresh_25.csv" ]; then

        if [ -f "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf" ]; then
        
            chrom_name=$(zgrep -v '^#' "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf" | awk 'NR==1 {print $1}')
        
            if [ "$chrom_name" = "Chromosome" ]; then
                bed_file="/home/sak0914/Mtb_Megapipe/references/WHO_catalog_resistance/regions.bed"
            else
                bed_file="/home/sak0914/Mtb_Megapipe/references/WHO_catalog_resistance/regions_refseq.bed"
            fi   
            
            bgzip -c "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf" > "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf.bgz"
                
            # tabix the bgzipped file, which will create fName.bgz.tbi
            tabix -p vcf -f "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf.bgz"
            
            bcftools view -R $bed_file "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf.bgz" | SnpSift extractFields '-' POS REF ALT FILTER QUAL IMPRECISE AF DP BQ MQ IC DC ANN -e "" > "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants.tsv"
            
            python3 -u /home/sak0914/Mtb_Megapipe/scripts/process_variants_for_WHO_catalog.py -i "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants.tsv"
            rm "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants.tsv"
            
            # get resistance predictions -- any Group 1 or 2 variant that passes QC leads to a prediction of R for a given drug. If not, predicted S
            python3 -u /home/sak0914/Mtb_Megapipe/scripts/WHO_catalog_resistance_pred.py -i "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_annot.tsv" -o "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_pred"
            python3 -u /home/sak0914/Mtb_Megapipe/scripts/WHO_catalog_resistance_pred.py -i "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_annot.tsv" -o "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_pred" --AF-thresh 0.25
    
            rm "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf.bgz"
            rm "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf.bgz.tbi"
    
            echo "Finished $sample_ID with $chrom_name"

        else
            echo "$out_dir/$sample_ID/WHO_resistance/${sample_ID}_variants_combinedCodons.eff.vcf doesn't exist"
        fi

    fi

done < "$input_file"