#! /usr/bin/env bash

nsnps=( "1k" "10k" "25k" "50k" "100k" "250k" "500k" )
for nsnp in "${nsnps[@]}"; do
    echo "./run_demuxalot_nsnp.sh ${nsnp}"
    echo "./run_demuxlet_nsnp.sh ${nsnp}"
    echo "./run_demux_vcf_nsnp.sh ${nsnp}"
    echo "./run_vireo_nsnp.sh ${nsnp}"
done

downsamp=( "0.1" "0.25" "0.5" "0.75" )
for ds in "${downsamp[@]}"; do
    echo "./run_demuxalot_downsample.sh ${ds}"
    echo "./run_demuxlet_downsample.sh ${ds}"
    echo "./run_demux_vcf_downsample.sh ${ds}"
    echo "./run_vireo_downsample.sh ${ds}"
done

ncell=( "100" "500" "1000" "5000" "10000" "20000" )
for nc in "${ncell[@]}"; do
    echo "./run_demuxalot_ncells.sh ${nc}"
    echo "./run_demuxlet_ncells.sh ${nc}"
    echo "./run_demux_vcf_ncells.sh ${nc}"
    echo "./run_vireo_ncells.sh ${nc}"
done

