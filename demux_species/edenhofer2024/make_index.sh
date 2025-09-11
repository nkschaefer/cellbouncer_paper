#! /usr/bin/env bash

# Get transcriptomes for human & macaque
wget https://ftp.ensembl.org/pub/release-114/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/macaca_fascicularis/cdna/Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa.gz

# Unpack
gunzip Homo_sapiens.GRCh38.cdna.all.fa.gz
gunzip Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa

# Prepare reference data
# Sample 20M most frequent unique k-mers from each transcriptome
${CELLBOUNCER}/utils/demux_species_ref.py -f Homo_sapiens.GRCh38.cdna.all.fa Macaca_fascicularis.Macaca_fascicularis_6.0.cdna.all.fa -n Human Macaque -N 20000000 -o HM20m

