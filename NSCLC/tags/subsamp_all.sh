#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: please provide path to CellRanger barcode whitelist (3M-february-2018.txt.gz)"
    >&2 echo "This is in [cellranger]/lib/python/cellranger/barcodes"
    exit 1
fi

bc=$1

# Make one set that is not downsampled
r1s=$( find . -mindepth 1 -maxdepth 2 -type f -name "*_R1_001.fastq.gz" | grep 40k_NSCLC_DTC_3p_HT_nextgem_cmo | sed 's/^/-1 /' | tr '\n' ' ' )
r2s=$( echo -e "$r1s" | sed 's/-1/-2/g' | sed 's/_R1_/_R2_/g' )
${CELLBOUNCER}/demux_tags $r1s $r2s -B filt_barcodes_translated.tsv -s seqs.txt -w $bc -o full

# Downsample reads

# Run all subsampling jobs
samps=( "0.001" "0.0025" "0.005" "0.01" "0.05" "0.1" "0.5" )
for dir_idx in $(seq 1 3 ); do
    dirn="40k_NSCLC_DTC_3p_HT_nextgem_cmo_${dir_idx}_fastqs"
    for samp in "${samps[@]}"; do
        for rep in $(seq 1 3 ); do
            ls ${dirn}/*_R1_001.fastq.gz | while read r1; do
                echo "./subsamp.sh ${r1} ${samp} ${rep}"
            done
        done
    done
done | parallel -j 12

# Use CellBouncer to count reads in all subsamples
find . -mindepth 1 -maxdepth 2 -type f -name "*_R1_001.fastq.gz" | grep -P "^s0" > r1s.txt

samps=( "0.001" "0.0025" "0.005" "0.01" "0.05" "0.1" "0.5" )
for samp in "${samps[@]}"; do
    for rep in $(seq 1 3 ); do
        name="s${samp}"
        if [ $rep -eq 2 ]; then
            name="${name}b"
        elif [ $rep -eq 3 ]; then
            name="${name}c"
        fi
        r1_this=$( cat r1s.txt | grep "$name" | sed 's/^/-1 /' | tr '\n' ' ' )
        r2_this=$( echo -e "$r1_this" | sed 's/-1/-2/g' | sed 's/_R1_/_R2_/g' )
        echo "${CELLBOUNCER}/demux_tags $r1_this $r2_this -B filt_barcodes_translated.tsv -s seqs.txt -w $bc -o $name"   
    done
done | parallel -j 12

# Replace barcodes with correct ones
for samp in "${samps[@]}"; do
    name="s${samp}"
    for rep in $(seq 1 3 ); do
        if [ $rep -eq 2 ]; then
            name="${name}b"
        elif [ $rep -eq 3 ]; then
            name="${name}c"
        fi
        ./swap_bcs.R "${name}.counts"
    done
done


