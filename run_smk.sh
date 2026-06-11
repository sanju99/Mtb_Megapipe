#!/bin/bash 
#SBATCH -c 8
#SBATCH -t 0-11:59
#SBATCH -p short
#SBATCH --mem=50G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu

source activate snakemake

snakemake --snakefile snakefile \
          --use-conda --conda-frontend conda --conda-prefix .snakemake/conda \
          --configfile config.yaml \
          --cores 8 --resources mem_mb=25000 \
          --rerun-incomplete --keep-going \
          --unlock


snakemake --snakefile snakefile \
          --use-conda --conda-frontend conda --conda-prefix .snakemake/conda \
          --configfile config.yaml \
          --cores 8 --resources mem_mb=25000 \
          --rerun-incomplete --keep-going
