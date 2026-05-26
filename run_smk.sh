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

snakemake --snakefile /home/sak0914/Mtb_Megapipe/snakefile \
          --use-conda --conda-frontend conda --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda \
          --configfile config.yaml \
          --directory /home/sak0914/Mtb_Megapipe \
          --cores 8 --resources mem_mb=25000 \
          --rerun-incomplete --keep-going \
          --unlock


snakemake --snakefile /home/sak0914/Mtb_Megapipe/snakefile \
          --use-conda --conda-frontend conda --conda-prefix /home/sak0914/Mtb_Megapipe/.snakemake/conda \
          --configfile config.yaml \
          --directory /home/sak0914/Mtb_Megapipe \
          --cores 8 --resources mem_mb=25000 \
          --rerun-incomplete --keep-going --allowed-rules pilon_get_variants_only_file exclude_low_confidence_sites combine_codon_variants annotate_variants_snpEff create_WHO_catalog_variants_TSV get_WHO_catalog_resistance_predictions #--dry-run