#!/bin/bash 
#SBATCH -c 1
#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

# source activate variant_calling_smk

# snakemake --snakefile snakefile_TRUST --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config_TRUST.yaml --cores 8 --directory /home/sak0914/MtbQuantCNN/analysis/TRUST --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --resources mem_mb=100000 --unlock

# snakemake --snakefile snakefile_TRUST --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config_TRUST.yaml --cores 8 --directory /home/sak0914/MtbQuantCNN/analysis/TRUST --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --resources mem_mb=100000

source activate bioinformatics

while IFS=$'\t', read -r sample_ID sample_ID_2 FQ_val
do

    if [ ! -f "/n/data1/hms/dbmi/farhat/rollingDB/TRUST/Illumina_culture_WGS_processed/$sample_ID/$sample_ID/kraken/kraken_report_standard_DB" ]; then

        echo $sample_ID

        kraken2 --db "/n/data1/hms/dbmi/farhat/mm774/References/Kraken2_DB_Dir/Kraken2_DB" --threads 8 --paired /n/data1/hms/dbmi/farhat/rollingDB/TRUST/Illumina_culture_WGS_processed/$sample_ID/$sample_ID/fastp/$sample_ID.R1.trimmed.fastq /n/data1/hms/dbmi/farhat/rollingDB/TRUST/Illumina_culture_WGS_processed/$sample_ID/$sample_ID/fastp/$sample_ID.R2.trimmed.fastq --report /n/data1/hms/dbmi/farhat/rollingDB/TRUST/Illumina_culture_WGS_processed/$sample_ID/$sample_ID/kraken/kraken_report_standard_DB > /n/data1/hms/dbmi/farhat/rollingDB/TRUST/Illumina_culture_WGS_processed/$sample_ID/$sample_ID/kraken/kraken_classifications_standard_DB

    fi

done < "/home/sak0914/Mtb_Megapipe/TRUST_samples.tsv"