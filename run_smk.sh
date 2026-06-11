#!/bin/bash 

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
