#! /usr/bin/env bash

props=( "0.05" "0.1" "0.15" "0.2" "0.25" )
seq 1 7 | while read d; do
    for p in ${props[@]}; do
        basen="d${d}c${p}"
        >&2 echo $basen
        cp ${basen}.counts ${basen}_corr.counts
        cp ${basen}.condf ${basen}_corr.condf
        cp ${basen}.samples ${basen}_corr.samples
        cp ../demux_vcf_sc_1000.assignments ${basen}_corr.assignments
        ${CELLBOUNCER}/quant_contam -o ${basen}_corr
    done
done
