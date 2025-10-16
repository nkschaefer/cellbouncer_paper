#! /usr/bin/env bash

for d in $(seq 1 7); do
    echo "../run_demuxalot.sh d${d}c0.05_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.05_demuxalot"
    echo "../run_demuxalot.sh d${d}c0.1_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.1_demuxalot"
    echo "../run_demuxalot.sh d${d}c0.15_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.15_demuxalot"
    echo "../run_demuxalot.sh d${d}c0.2_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.2_demuxalot"
    echo "../run_demuxalot.sh d${d}c0.25_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.25_demuxalot"
done
