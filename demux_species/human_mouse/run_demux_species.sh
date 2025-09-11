#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: provide k"
    exit 1
fi

k=$1

if [ -d kmer_${k}_cellbouncer ]; then
    rm -r kmer_${k}_cellbouncer
fi

mprof run -o mprofile_run_kmer_${k}.dat ${CELLBOUNCER}/demux_species -k hm${k} -r 5k_hgmm_3p_nextgem_fastqs/ALL_R1.fastq.gz -R 5k_hgmm_3p_nextgem_fastqs/ALL_R2.fastq.gz -d -o kmer_${k}_cellbouncer -w wl.txt 
