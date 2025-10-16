#! /usr/bin/env bash

for d in $(seq 1 7); do
    echo "../run_vireo.sh d${d}c0.05_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.05_vireo"
    echo "../run_vireo.sh d${d}c0.1_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.1_vireo"
    echo "../run_vireo.sh d${d}c0.15_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.15_vireo"
    echo "../run_vireo.sh d${d}c0.2_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.2_vireo"
    echo "../run_vireo.sh d${d}c0.25_merged.bam ../samp1000.barcodes2 ../500k.vcf.gz d${d}c0.25_vireo"
done
