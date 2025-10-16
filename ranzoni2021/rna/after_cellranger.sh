#! /usr/bin/env bash
if [ $# -lt 1 ]; then
    >&2 echo "ERROR: please provide cellranger output dir"
    exit 1
fi

# First, refine VCF
${CELLBOUNCER}/utils/refine_vcf -b ${1}/outs/possorted_genome_bam.bam -v vars.filt.vcf.gz -p 0.01 -a cellranger_bc2id.assignments > vars.filt.refined.vcf.gz
tabix -p vcf vars.filt.refined.vcf.gz

# Now run demux_vcf
${CELLBOUNCER}/demux_vcf -v vars.filt.refined.vcf.gz -b ${1}/outs/possorted_genome_bam.bam -o demux_vcf_withdoub
${CELLBOUNCER}/demux_vcf -v vars.filt.refined.vcf.gz -b ${1}/outs/possorted_genome_bam.bam -D 0 -o demux_vcf_nodoub

# Run quant_contam
${CELLBOUNCER}/quant_contam -o demux_vcf_withdoub
${CELLBOUNCER}/quant_contam -o demux_vcf_nodoub -D 0 

# Run DecontX
./decontX.R $1

# Plot comparison
./make_comp_plot.R

