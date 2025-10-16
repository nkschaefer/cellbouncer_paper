#! /usr/bin/env bash

for d in $(seq 1 7); do
    echo "../run_demuxlet.sh d${d}c0.05_merged.bam ../500k.vcf.gz d${d}c0.05_demuxlet"
    echo "../run_demuxlet.sh d${d}c0.1_merged.bam ../500k.vcf.gz d${d}c0.1_demuxlet"
    echo "../run_demuxlet.sh d${d}c0.15_merged.bam ../500k.vcf.gz d${d}c0.15_demuxlet"
    echo "../run_demuxlet.sh d${d}c0.2_merged.bam ../500k.vcf.gz d${d}c0.2_demuxlet"
    echo "../run_demuxlet.sh d${d}c0.25_merged.bam ../500k.vcf.gz d${d}c0.25_demuxlet"
done
