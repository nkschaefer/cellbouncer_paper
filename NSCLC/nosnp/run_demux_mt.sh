#! /usr/bin/env bash

# Run CellBouncer
${CELLBOUNCER}/demux_mt -c -b ../40k_all.bam -o demux_mt_all
${CELLBOUNCER}/demux_mt -c -b ../40k_sc_2000.plusjunk.chosen.bam -o demux_mt_2000
${CELLBOUNCER}/demux_mt -c -b ../40k_sc_5000.plusjunk.chosen.bam -o demux_mt_5000

# Make plot
${CELLBOUNCER}/plot/demux_mt_clust demux_mt_all
${CELLBOUNCER}/plot/demux_mt_clust demux_mt_2000
${CELLBOUNCER}/plot/demux_mt_clust demux_mt_5000

# Insert read groups
${CELLBOUNCER}/utils/bam_indv_rg -b ../40k_all.bam -a demux_mt_all.assignments -o 40k_all_rg.bam
${CELLBOUNCER}/utils/bam_indv_rg -b ../40k_sc_2000.plusjunk.chosen.bam -a demux_mt_2000.assignments -o 40k_sc_2000_rg.bam
${CELLBOUNCER}/utils/bam_indv_rg -b ../40k_sc_5000.plusjunk.chosen.bam -a demux_mt_5000.assignments -o 40k_sc_5000_rg.bam
samtools index 40k_all_rg.bam
samtools index 40k_sc_2000_rg.bam
samtools index 40k_sc_5000_rg.bam

# Run FreeBayes
cat ../hg38.fa.fai | cut -f1-2 | grep -v chrUn | grep -v random | grep -v chrM | grep -v chrY | grep -v alt | while read line; do
    chrom=$( echo -e "$line" | cut -f1 )
    chromlen=$( echo -e "$line" | cut -f2 )
    echo "freebayes -f ../hg38.fa -L 40k_all_rg.bam -r ${chrom}:0-${chromlen} > vars.all.${chrom}.vcf" 
    echo "freebayes -f ../hg38.fa -L 40k_sc_2000_rg.bam -r ${chrom}:0-${chromlen} > vars.2000.${chrom}.vcf"
    echo "freebayes -f ../gh38.fa -L 40k_sc_5000_rg.bam -r ${chrom}:0-${chromlen} > vars.5000.${chrom}.vcf"
done | parallel -j 12

cat vars.all.chr21.vcf | grep -P "^#" > MT_all.vcf
cat vars.2000.chr21.vcf | grep -P "^#" > MT_2000.vcf
cat vars.5000.chr21.vcf | grep -P "^#" > MT_5000.vcf

seq 1 22 | while read chrom; do
    cat vars.all.chr${chrom}.vcf | grep -v -P "^#" >> MT_all.vcf
    cat vars.2000.chr${chrom}.vcf | grep -v -P "^#" >> MT_2000.vcf
    cat vars.5000.chr${chrom}.vcf | grep -v -P "^#" >> MT_5000.vcf
done
cat vars.all.chrX.vcf | grep -v -P "^#" >> MT_all.vcf
cat vars.2000.chrX.vcf | grep -v -P "^#" >> MT_2000.vcf
cat vars.5000.chrX.vcf | grep -v -P "^#" >> MT_5000.vcf

bgzip MT_all.vcf
bgzip MT_2000.vcf
bgzip MT_5000.vcf
tabix -p vcf MT_all.vcf
tabix -p vcf MT_2000.vcf
tabix -p vcf MT_5000.vcf

# Filter variants
bcftools view -c 1 -C 13 -v snps -i 'QUAL > 100 & F_MISSING < 0.5' -O z MT_all.vcf.gz > MT_all.filt.vcf.gz
# Use more lenient filtering for downsampled data
bcftools view -c 1 -C 13 -v snps -i 'QUAL > 50 & F_MISSING < 0.75' -O z MT_2000.vcf.gz > MT_2000.filt.vcf.gz
bcftools view -c 1 -C 13 -v snps -i 'QUAL > 50 & F_MISSING < 0.75' -O z MT_5000.vcf.gz > MT_5000.filt.vcf.gz
tabix -p vcf MT_all.filt.vcf.gz
tabix -p vcf MT_2000.filt.vcf.gz
tabix -p vcf MT_5000.filt.vcf.gz

# Run demux_vcf
${CELLBOUNCER}/demux_vcf -v MT_all.filt.vcf.gz -b ../40k_all.bam -o MTvars_all
${CELLBOUNCER}/demux_vcf -v MT_2000.filt.vcf.gz -b ../40k_sc_2000.plusjunk.chosen.bam -o mitovars_2000
${CELLBOUNCER}/demux_vcf -v MT_5000.filt.vcf.gz -b ../40k_sc_5000.plusjunk.chosen.bam -o mitovars_5000

# Refine genotypes
${CELLBOUNCER}/utils/refine_vcf -b ../40k_all.bam -v MT_all.filt.vcf.gz -a MTvars_all.assignments -p 0 > MTvars_all_refined.vcf.gz
${CELLBOUNCER}/utils/refine_vcf -b ../40k_sc_2000.plusjunk.chosen.bam -v MT_2000.filt.vcf.gz -a mitovars_2000.assignments -p 0 > mitovars_2000_refined.vcf.gz
${CELLBOUNCER}/utils/refine_vcf -b ../40k_sc_5000.plusjunk.chosen.bam -v MT_5000.filt.vcf.gz -a mitovars_5000.assignments -p 0 > mitovars_5000_refined.vcf.gz
tabix -p vcf MTvars_all_refined.vcf.gz
tabix -p vcf mitovars_2000_refined.vcf.gz
tabix -p vcf mitovars_5000_refined.vcf.gz

# Run demux_vcf again
${CELLBOUNCER}/demux_vcf -v MTvars_all_refined.vcf.gz -b ../40k_all.bam -o MTvars_all_refined
${CELLBOUNCER}/demux_vcf -v mitovars_2000_refined.vcf.gz -b ../40k_sc_2000.plusjunk.chosen.bam -o mitovars_2000_refined
${CELLBOUNCER}/demux_vcf -v mitovars_5000_refined.vcf.gz -b ../40k_sc_5000.plusjunk.chosen.bam -o mitovars_5000_refined





