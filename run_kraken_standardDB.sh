#!/bin/bash 
#SBATCH -c 4
#SBATCH -t 0-11:59
#SBATCH -p short 
#SBATCH --mem=50G
#SBATCH -o /home/sak0914/Errors/zerrors_%j.out 
#SBATCH -e /home/sak0914/Errors/zerrors_%j.err 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=skulkarni@g.harvard.edu


#################### REQUIRED COMMAND LINE ARGUMENTS (IN THIS ORDER): ####################

# 1. TSV file with 2 columns (AND NO HEADER): 
    # col1: sample_ID, i.e. "SAMEA104362063"
    # col2: comma-separated string of all sequencing runs for the sample in col1, i.e. "ERR2184222,ERR9029914"
# 2. fastq_dir: directory where paired-end Illumina sequencing reads are stored
# 3. out_dir: output directory where results should be stored. i.e. /n/data1/hms/dbmi/farhat/rollingDB/cryptic_output or /n/data1/hms/dbmi/farhat/rollingDB/genomic_data

##### NOTE: The values in col1 and col2 of the input file can be the same, i.e. for unpublished data. But for standardization, it is best to use the BioSample accession ID in col1 and sequencing run IDs in col2 #####

set -o errexit # any error will cause the shell script to exit immediately. This is not native bash behavior
source activate bioinformatics # CHANGE TO YOUR OWN ENVIRONMENT OR REMOVE IF RUNNING IN THE BASE ENV

if ! [ $# -eq 2 ]; then
    echo "Please pass in 2 command line arguments: a text file with two columns: sample_ID and sequencing runs, and the output directory"
    exit
fi

# command line arguments
input_file=$1
out_dir=$2 # directory where completed variant calling results for sample_ID IDs will be stored (OUTPUT)

# Read the TSV file line by line, skiping the header. IFS sets the field separator. Here, it is tab
while IFS=$'\t', read -r sample_ID run_IDs_string download_FQ
do

    if [[ -f "$out_dir/$sample_ID/$sample_ID/fastp/$sample_ID.R1.trimmed.fastq" && -f "$out_dir/$sample_ID/$sample_ID/fastp/$sample_ID.R2.trimmed.fastq" ]]; then

        echo "Performing Kraken classification with the standard database on $sample_ID"

        kraken2 --db "/n/data1/hms/dbmi/farhat/mm774/References/Kraken2_DB_Dir/Kraken2_DB" --threads 8 --paired "$out_dir/$sample_ID/$sample_ID/fastp/$sample_ID.R1.trimmed.fastq" "$out_dir/$sample_ID/$sample_ID/fastp/$sample_ID.R2.trimmed.fastq" --report "$out_dir/$sample_ID/$sample_ID/kraken/kraken_report_standard_DB" > "$out_dir/$sample_ID/$sample_ID/kraken/kraken_classifications_standard_DB"

    else

        echo "No fastp-trimmed FASTQ files yet for $sample_ID"
        
    fi

done < "$input_file"