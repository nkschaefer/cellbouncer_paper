#! /usr/bin/env bash

wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/068/E-MTAB-9068/Files/E-MTAB-9068.sdrf.txt
./download_seqs.sh

# Count cells per sample
ncell1=$( cat E-MTAB-9068.sdrf.txt | grep -w Sample1 | grep "_1.fastq.gz" | wc -l )
ncell2=$( cat E-MTAB-9068.sdrf.txt | grep -w Sample2 | grep "_1.fastq.gz" | wc -l )
ncell3=$( cat E-MTAB-9068.sdrf.txt | grep -w Sample3 | grep "_1.fastq.gz" | wc -l )

mkdir Sample1_bc
mkdir Sample2_bc
mkdir Sample3_bc

# Add barcodes
seq 1 $ncell1 | while read idx; do
    echo "./addbc_sample.sh Sample1 $idx"
done | parallel -j 12
seq 1 $ncell2 | while read idx; do
    echo "./addbc_sample.sh Sample2 $idx"
done | parallel -j 12
seq 1 $ncell3 | while read idx; do
    echo "./addbc_sample.sh Sample3 $idx"
done | parallel -j 12

# Get reference genome & minimap2 idx
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
minimap2 -d hg38.idx hg38.fa

# Align data
./map_sample1.sh 
./map_sample2.sh
./map_sample3.sh

# Extract chrM from each
samtools view -bh Sample1_sorted.bam chrM > Sample1_chrM.bam
samtools view -bh Sample2_sorted.bam chrM > Sample2_chrM.bam
samtools view -bh Sample3_sorted.bam chrM > Sample3_chrM.bam

samtools index Sample1_chrM.bam
samtools index Sample2_chrM.bam
samtools index Sample3_chrM.bam

# Combine into one chrM file
samtools merge chrM.bam Sample1_chrM.bam Sample2_chrM.bam Sample3_chrM.bam
samtools index chrM.bam

# Infer mitochondrial haplotypes
${CELLBOUNCER}/demux_mt -b chrM.bam -o chrM_demux_mt

# Plot
${CELLBOUNCER}/plot/demux_mt_clust.R chrM_demux_mt




