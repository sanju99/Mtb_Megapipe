import numpy as np
import pandas as pd
import glob, os, argparse, time
from Bio import Seq, SeqIO, Entrez

parser = argparse.ArgumentParser()

# "/n/data1/hms/dbmi/farhat/rollingDB/TRUST/clinical_data/20240826_raw_data.csv"
parser.add_argument("-t", "--taxid", dest='taxid', type=int, required=True, help='Taxid for which to get reads')
parser.add_argument("-d", "--database", dest='db_name', type=str, required=True, help='Database name from which to get taxonomy')
parser.add_argument("-i", "--input", dest="in_fName", type=str, required=True, help='Kraken classifications file')

cmd_line_args = parser.parse_args()
taxid = cmd_line_args.taxid
db_name = cmd_line_args.db_name
in_fName = cmd_line_args.in_fName
out_dir = os.path.dirname(in_fName)


def load_nodes_dmp(nodes_path):
    parent_map = {}
    child_map = {}
    with open(nodes_path, 'r') as f:
        for line in f:
            parts = [x.strip() for x in line.split('|')]
            child_taxid, parent_taxid = int(parts[0]), int(parts[1])
            parent_map[child_taxid] = parent_taxid
            if parent_taxid not in child_map:
                child_map[parent_taxid] = []
            child_map[parent_taxid].append(child_taxid)
    return parent_map, child_map
    

def get_parent_taxids(taxid, parent_map):
    ancestors = []
    while taxid != 1:  # 1 is the root (cellular organisms)
        taxid = parent_map[taxid]
        ancestors.append(taxid)
    return ancestors


def get_child_taxids(taxid, child_map):
    descendants = []
    stack = [taxid]
    while stack:
        current = stack.pop()
        children = child_map.get(current, [])
        descendants.extend(children)
        stack.extend(children)
    return descendants


# analyze taxonomy map 
parent_map, child_map = load_nodes_dmp(f"{db_name}/taxonomy/nodes.dmp")

# keep only the child taxids
child_taxids = get_child_taxids(taxid, child_map)

blastn_columns = [
    "qseqid",    # Query sequence ID
    "sseqid",    # Subject (database) sequence ID
    "pident",    # Percentage of identical matches
    "length",    # Alignment length
    "mismatch",  # Number of mismatches
    "gapopen",   # Number of gap openings
    "qstart",    # Start of alignment in query
    "qend",      # End of alignment in query
    "sstart",    # Start of alignment in subject
    "send",      # End of alignment in subject
    "evalue",    # Expectation value (E-value)
    "bitscore"   # Bit score
]

df_blast_results = pd.read_csv(in_fName, sep="\t", header=None, names=blastn_columns).sort_values("bitscore", ascending=False)

# separate the taxid, which is at the end of the query genome name (because that's how we made it for the kraken database)
df_blast_results['TaxID'] = df_blast_results['sseqid'].str.split('|').str[-1].astype(int)

# get the top taxon per read
top_taxon_per_read = df_blast_results.drop_duplicates('qseqid')[['qseqid', 'TaxID']]

# extract and save
keep_reads = top_taxon_per_read.query("TaxID in @child_taxids")['qseqid']

print(f"Keeping {len(keep_reads)}/{len(top_taxon_per_read)} ({np.round(len(keep_reads)/len(top_taxon_per_read)*100)}%) reads")

keep_reads.to_csv(f"{out_dir}/MTBC_read_names.txt", sep='\t', header=None, index=False)