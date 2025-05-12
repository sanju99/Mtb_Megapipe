

# determine the appropriate FASTA file by getting the name of the reference chromosome used for alignment
chrom_name=$(samtools idxstats $bam_file  | cut -f1 | head -1)

if [ "$chrom_name" = "Chromosome" ]; then
    ref_fasta="/home/sak0914/Mtb_Megapipe/references/ref_genome/H37Rv_NC_000962.3.fna"
else
    ref_fasta="/home/sak0914/Mtb_Megapipe/references/ref_genome/refseq.fna"
fi

bam_file="SAMEA104362043.dedup.bam"
cram_file="SAMEA104362043.dedup.cram"

samtools view -T $ref_fasta -C -o $cram_file $bam_file

samtools view -T $ref_fasta -b -o $new_bam_file $cram_file

samtools view $new_bam_file > $new_bam_txt
samtools view $cram_file > $cram_txt

if cmp -s bam2.txt cram.txt; then
    echo "Files are identical"
else
    echo "Files differ"
    exit 1
fi