#!/bin/bash 
#SBATCH -c 1
#SBATCH -t 1-00:00
#SBATCH -p medium
#SBATCH --mem=10G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

source activate snakemake

snakemake --snakefile /home/sak0914/Mtb_Megapipe/snakefile_cleaning --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config_2.yaml --directory /home/sak0914/who-analysis --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --cores 8 --resources mem_mb=5000 --unlock

snakemake --snakefile /home/sak0914/Mtb_Megapipe/snakefile_cleaning --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config_2.yaml  --directory /home/sak0914/who-analysis --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --cores 8 --resources mem_mb=5000 #--dry-run #--rerun-triggers mtime