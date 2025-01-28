#!/bin/bash 
#SBATCH -c 4
#SBATCH -t 4-23:59
#SBATCH -p medium
#SBATCH --mem=50G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

source activate variant_calling_smk

snakemake --snakefile snakefile --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config.yaml --cores 8 --directory /home/sak0914/Mtb_Megapipe --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --resources mem_mb=10000 --unlock

snakemake --snakefile snakefile --use-conda --conda-frontend conda --rerun-incomplete --keep-going --configfile config.yaml --cores 8 --directory /home/sak0914/Mtb_Megapipe --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda --resources mem_mb=10000