#! /usr/bin/env bash

# Download and unpack Ensembl annotations
wget https://ftp.ensembl.org/pub/release-114/fasta/macaca_mulatta/cdna/Macaca_mulatta.Mmul_10.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/tupaia_belangeri/cdna/Tupaia_belangeri.TREESHREW.cdna.all.fa.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/rattus_norvegicus/cdna/Rattus_norvegicus.GRCr8.cdna.all.fa.gz

gunzip Macaca_mulatta.Mmul_10.cdna.all.fa.gz
gunzip Tupaia_belangeri.TREESHREW.cdna.all.fa.gz
gunzip Rattus_norvegicus.GRCr8.cdna.all.fa.gz

# Download all SRA data
cat sra_samples.tsv | tail -n +2 | grep SRR | cut -f2 | sort | uniq | while read samp; do
    if [ ! -d $samp ]; then
        mkdir $samp
    fi
done

# Download samples in parallel (uses 8 cores)
# Note: there seems to be a limit of 10 simultaneous SRA downloads
cat sra_samples.tsv | tail -n +2 | grep SRR | while read line; do
    run=$( echo -e "$line" | cut -f1 )
    samp=$( echo -e "$line" | cut -f2 )
    echo "cd $samp && fastq-dump --split-3 $run"
done | parallel -j 8

# GZIP all FASTQs
find . -mindepth 2 -maxdepth 2 -type f -name "*.fastq" | grep SRR | while read fn; do
    echo "gzip ${fn}" 
done | parallel -j 8

# Set up symlinks for demux_species pipeline
./make_symlinks.sh

