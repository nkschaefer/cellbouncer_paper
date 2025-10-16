#! /usr/bin/env bash

props=( "0.05" "0.1" "0.15" "0.2" "0.25" )
seq 1 7 | while read d; do
    for p in ${props[@]}; do
        basen="d${d}c${p}"
        >&2 echo $basen
        ${CELLBOUNCER}/bulkprops -e 0.001 -N 0 -v ../500k.vcf.gz -o $basen -b ${basen}_merged.bam 
    done
done
