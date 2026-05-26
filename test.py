import numpy as np
import pandas as pd
import os, vcf, sys, itertools, shutil, subprocess, glob, yaml
from Bio import Entrez, Seq, SeqIO
import warnings
warnings.filterwarnings("ignore")
import seaborn as sns
import matplotlib.pyplot as plt

H37Rv_regions = pd.read_csv("./references/ref_genome/mycobrowser_h37rv_v4.csv")
isolate_metadata = pd.read_csv("/n/data1/hms/dbmi/farhat/rollingDB/metadata/isolate_metadata.csv")
genomic_data_dir = "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data"

def extract_kraken_reports(df, sample_id_col, run_id_col, out_dir=None):

    df_add = pd.DataFrame(columns = [sample_id_col, run_id_col, 'Kraken_Unclassified_Percent'])
    k = 0

    if out_dir is None:
        out_dir = "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data"
    
    for i, name in enumerate(df[sample_id_col].values):

        for run_id in np.sort(df.query(f"{sample_id_col}==@name")[run_id_col].values):
    
            if os.path.isfile(f"{out_dir}/{name}/{run_id}/kraken/kraken_report"):
        
                df_kraken = pd.read_csv(f"{out_dir}/{name}/{run_id}/kraken/kraken_report", sep='\t', header=None)

                kraken_unclassified = df_kraken.loc[df_kraken[3]=='U'][0].values[0]

                df_add.loc[k, :] = [name, run_id, kraken_unclassified]
                k += 1
    
    return df_add
    
    
    
    
def compute_BAM_depth_metrics(df, sample_id_col, run_id_col, out_dir=None):

    df_BAM_depths = pd.DataFrame(columns = [sample_id_col, run_id_col, 'Mean_Depth', 'Median_Depth', 'Prop_20x', 'Prop_10x'])
    idx = 0

    if out_dir is None:
        out_dir = "/n/data1/hms/dbmi/farhat/rollingDB/genomic_data"
    
    for i, name in enumerate(df[sample_id_col].values):

        run_ids = np.sort(df.query(f"{sample_id_col}==@name")[run_id_col].values)
    
        if os.path.isfile(f"{out_dir}/{name}/bam/{name}.depth.tsv.gz") and len(run_ids) > 0:
    
            df_depth = pd.read_csv(f"{out_dir}/{name}/bam/{name}.depth.tsv.gz", compression='gzip', header=None, sep='\t')

            pass_props = []
            
            # the first 2 columns are CHROM and POS. The remaining columns are the depth for each run
            # assert len(run_ids) == df_depth.shape[1]-2

            for k, col in enumerate(df_depth.columns[2:]):

                mean_depth = df_depth[col].mean()
                median_depth = df_depth[col].median()
                prop_20x = len(df_depth.loc[df_depth[col] >= 20]) / len(df_depth)
                prop_10x = len(df_depth.loc[df_depth[col] >= 10]) / len(df_depth)

                try:
                    df_BAM_depths.loc[idx, :] = [name, run_ids[k], mean_depth, median_depth, prop_20x, prop_10x]
                    idx += 1    
                except:
                    print(f"{name} failed")
        else:
            print(f"No depth file for {name}")
    
        if i % 1000 == 0:
            pd.concat([df_BAM_finished, df_BAM_depths]).to_csv(out_file, index=False)
            print(i)
    
    return df_BAM_depths


out_file = "Anna_samples_BAM_depth_summary.csv"

if os.path.isfile(out_file):
    df_BAM_finished = pd.read_csv(out_file)
    print(f"Already finished checking {df_BAM_finished.ROLLINGDB_ID.nunique()} files")
else:
    df_BAM_finished = pd.DataFrame(columns = ['ROLLINGDB_ID'])

df_Anna = pd.read_csv("~/data_cleaning/data/Anna_binary_resistance_table_to_rerun.csv")

df_Anna_metadata = isolate_metadata.query("ROLLINGDB_ID in @df_Anna.New_ID & ROLLINGDB_ID not in @df_BAM_finished.ROLLINGDB_ID")[['ROLLINGDB_ID', 'RUN']].drop_duplicates()
print(f"{len(df_Anna_metadata)} files left to check")

df_BAM_depths = compute_BAM_depth_metrics(df_Anna_metadata, "ROLLINGDB_ID", 'RUN', out_dir=None)

pd.concat([df_BAM_finished, df_BAM_depths]).drop_duplicates().to_csv(out_file, index=False)