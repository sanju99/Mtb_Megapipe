import numpy as np
import pandas as pd
import glob, os, argparse, time
from Bio import Seq, SeqIO, Entrez

parser = argparse.ArgumentParser()

# "/n/data1/hms/dbmi/farhat/rollingDB/TRUST/clinical_data/20240826_raw_data.csv"
parser.add_argument("-t", "--taxid", dest='taxid', type=int, required=True, help='Taxid for which to get reads')
parser.add_argument("-d", "--database", dest='db_name', type=str, required=True, help='Database name from which to get taxonomy')
# parser.add_argument("-e", "--email", dest='email', type=str, required=True, help='Email address for querying Entrez tools')
parser.add_argument("-i", "--input", dest="in_fName", type=str, required=True, help='Kraken classifications file')
# parser.add_argument("-o", "--output", dest='out_fName', type=str, required=True, help='File of read names to keep')
parser.add_argument("--include-children", dest='include_children', action='store_true', help='Include children taxa of the given taxid')
parser.add_argument("--include-parents", dest='include_parents', action='store_true', help='Include parent taxa of the given taxid')

cmd_line_args = parser.parse_args()
taxid = cmd_line_args.taxid
db_name = cmd_line_args.db_name
# email = cmd_line_args.email
include_children = cmd_line_args.include_children
include_parents = cmd_line_args.include_parents

in_fName = cmd_line_args.in_fName
out_dir = os.path.dirname(in_fName)
# out_fName = cmd_line_args.out_fName


# def get_child_taxids(parent_taxid, email, max_results=500):
#     '''
#     1000 is the largest number of requests that can be returned at once. 
    
#     Here, we're going to chunk it up into groups of 500 so that we don't submit too many requests at once.
#     '''
    
#     Entrez.email = email

#     # Initialize variables. Include the parent taxid
#     all_taxids = [parent_taxid]
#     ret_start = 0

#     while True:
        
#         # Search for child taxids and keep going until all have been saved to a list
#         search_handle = Entrez.esearch(db="taxonomy", term=f"txid{parent_taxid}[Subtree]", retstart=ret_start, retmax=max_results)
#         search_results = Entrez.read(search_handle)
#         search_handle.close()

#         # Extract the list of taxids from the current batch
#         id_list = search_results["IdList"]
        
#         if not id_list:
#             break  # No more results to fetch

#         # Fetch details for the found taxids
#         fetch_handle = Entrez.efetch(db="taxonomy", id=",".join(id_list), retmode="xml")
#         fetch_results = Entrez.read(fetch_handle)
#         fetch_handle.close()

#         # Extract the taxids and add to cumulative list
#         child_taxids = [taxon["TaxId"] for taxon in fetch_results]
#         all_taxids.extend(np.unique(child_taxids))

#         # if the number of returned taxids is the same as the maximum that can be returned, it means there are more than can be retrieved
#         if len(child_taxids) == max_results:

#             # so increment ret_start and keep going
#             ret_start += max_results

#             # print progress
#             print(f"Found {max_results} results")

#             # wait to prevent a too many requests error
#             time.sleep(1)

#         # break out of the loop if we've retrieved all results
#         else:
#             break

#     return np.array(all_taxids).astype(int)




# def get_parent_taxids(taxid, email):
#     """
#     Given a taxid, return all its parent taxids (up to the root of the NCBI taxonomy).
#     Includes the input taxid itself as the first element.
#     """
#     Entrez.email = email

#     # Fetch taxonomy information for the given taxid
#     handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
#     records = Entrez.read(handle)
#     handle.close()

#     if not records:
#         raise ValueError(f"No taxonomy record found for taxid {taxid}")

#     lineage = records[0].get("LineageEx", [])
#     parent_taxids = [int(taxid)]  # Include the input taxid

#     # Add the taxids from the lineage
#     parent_taxids.extend(int(ancestor["TaxId"]) for ancestor in lineage)

#     return np.array(parent_taxids, dtype=int)


# # include the argument taxid
# select_taxids = [taxid]

# if include_children:
#     child_taxids = list(get_child_taxids(taxid, email))
# else:
#     child_taxids = []

# if include_parents:
#     parent_taxids = list(get_parent_taxids(taxid, email))
# else:
#     parent_taxids = []

# # combine into a single list
# select_taxids += child_taxids + parent_taxids



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

# # use the functions to get all taxids
# if include_parents:
#     parent_taxids = get_parent_taxids(taxid, parent_map)
# else:
#     parent_taxids = []

# if include_children:
#     child_taxids = get_child_taxids(taxid, child_map)
# else:
#     child_taxids = []

# use the functionss to get all taxids
parent_taxids = get_parent_taxids(taxid, parent_map)
child_taxids = get_child_taxids(taxid, child_map)

df_kraken_classifications = pd.read_csv(in_fName, sep='\t', header=None)
df_kraken_classifications.columns = ['Classified', 'ReadName', 'TaxID', 'Length', 'LCA']

# keep the argument taxid and children
select_taxids = [taxid] + child_taxids
keep_reads = df_kraken_classifications.query("TaxID in @select_taxids")['ReadName']

# blast the parent taxids
blast_reads = df_kraken_classifications.query("TaxID in @parent_taxids")['ReadName']

print(f"Keeping {len(keep_reads)}/{len(df_kraken_classifications)} ({np.round(len(keep_reads)/len(df_kraken_classifications)*100)}%) reads")
print(f"Writing {len(blast_reads)}/{len(df_kraken_classifications)} ({np.round(len(blast_reads)/len(df_kraken_classifications)*100)}%) reads to FASTA for BLAST")

keep_reads.to_csv(f"{out_dir}/MTBC_read_names.txt", sep='\t', header=None, index=False)
blast_reads.to_csv(f"{out_dir}/BLAST_read_names.txt", sep='\t', header=None, index=False)