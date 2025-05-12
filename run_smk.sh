#!/bin/bash 
#SBATCH -c 1
#SBATCH -t 2-23:59
#SBATCH -p medium
#SBATCH --mem=10G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

source activate snakemake

snakemake --snakefile snakefile --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config_2.yaml --config isolates_to_run="/home/sak0914/Mtb_Megapipe/BAM_files_conversion.tsv" --directory /home/sak0914/who-analysis --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --cores 8 --resources mem_mb=5000 --unlock

snakemake --snakefile snakefile --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config_2.yaml --config isolates_to_run="/home/sak0914/Mtb_Megapipe/BAM_files_conversion.tsv" --directory /home/sak0914/who-analysis --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --cores 8 --resources mem_mb=5000 #--rerun-triggers mtime