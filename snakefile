import numpy as np
import pandas as pd

configfile: "config.yaml"

# Define the list of samples, either directly or through config
output_dir = config["output_dir"]

# dataframe of isolates to run. First column is the sample ID, second column is a comma-separated string of the sequencing run IDs
df_samples_runs = pd.read_csv(config['isolates_to_run'], sep='\t', header=None)
sample_run_dict = dict(zip(df_samples_runs[0], df_samples_runs[1].str.split(',')))

# make a dictionary of binary values indicating whether (1) or not (0) to download publicly available FASTQ files
assert df_samples_runs[2].isin([0, 1]).all()
download_public_FASTQ_dict = dict(zip(df_samples_runs[0], df_samples_runs[2]))

include: "rules/cleaning.smk"

# Define rules in the Snakefile

rule all:
    input:
        [f"{output_dir}/{sample_ID}/bam/{sample_ID}.dedup.cram" for sample_ID in sample_run_dict.keys()],
        # [f"{output_dir}/{sample_ID}/{run_ID}/bam/{run_ID}.dedup.cram" for sample_ID in sample_run_dict.keys() for run_ID in sample_run_dict[sample_ID]],