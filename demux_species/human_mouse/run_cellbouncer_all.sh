#! /usr/bin/env bash

# Download transcriptomes (from GENCODE)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.transcripts.fa.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M38/gencode.vM38.transcripts.fa.gz

# Unpack
gunzip gencode.v49.transcripts.fa.gz
gunzip gencode.vM38.transcripts.fa.gz

./run_index.sh 22
./run_demux_species.sh 22
./run_index.sh 32
./run_demux_species.sh 32
./run_index.sh 42 
./run_demux_species.sh 42
./run_index.sh 52
./run_demux_species.sh 52


