#! /usr/bin/env bash

# Run bulkprops on full data set
${CELLBOUNCER}/bulkprops -b ../40k_all.bam -v ../500k.vcf -e 0.001 -o demux_vcf_donor_all_500k

# Run once for each program independently
${CELLBOUNCER}/doublet_dragon demux_vcf_dd ../demux_vcf_donor_all_500k.assignments
${CELLBOUNCER}/doublet_dragon demux_tags_dd ../tags/full.assignments
${CELLBOUNCER}/doublet_dragon demux_mt_dd ../nosnp/demux_mt_all.assignments

# Create the plot
./make_plot.R

