#! /usr/bin/env bash

# Create a set of "naive" variants for use with programs that
# require prior knowledge of variant sites (call as a pool)

freebayes -f ../hg38.fa -b ../40k_all.bam > varsNaive.vcf
bgzip varsNaive.vcf
tabix -p vcf varsNaive.vcf.gz

freebayes -f ../hg38.fa -b ../40k_sc_2000.chosen.bam > varsNaive2000.vcf
bgzip varsNaive2000.vcf
tabix -p vcf varsNaive2000.vcf.gz

freebayes -f ../hg38.fa -b ../40k_sc_5000.chosen.bam > varsNaive5000.vcf
bgzip varsNaive5000.vcf
tabix -p vcf varsNaive5000.vcf.gz

bcftools view -v snps -m 2 -M 2 -i 'QUAL > 100 & F_MISSING < 0.5' -O z varsNaive.vcf.gz > varsNaive.filt.q100.vcf.gz
tabix -p vcf varsNaive.filt.q100.vcf.gz

bcftools view -v snps -m 2 -M 2 -i 'QUAL > 50 & F_MISSING < 0.75' -O z varsNaive2000.vcf.gz > varsNaive2000.filt.q50.vcf.gz
tabix -p vcf varsNaive2000.filt.q50.vcf.gz

bcftools view -v snps -m 2 -M 2 -i 'QUAL > 50 & F_MISSING < 0.75' -O z varsNaive5000.vcf.gz > varsNaive5000.filt.q50.vcf.gz
tabix -p vcf varsNaive5000.filt.q50.vcf.gz

# Run all Cellbouncer stuff
./run_demux_mt.sh

# Run freemuxlet
./run_freemuxlet.sh ../40k_all.bam ../donor_all.bcs varsNaive.filt.q100.vcf.gz freemuxlet_all
./run_freemuxlet.sh ../40k_sc_2000.chosen.bam ../samp2000.barcodes2 varsNaive2000.filt.q50.vcf.gz freemuxlet_2000
./run_freemuxlet.sh ../40k_sc_5000.chosen.bam ../samp5000.barcodes2 varsNaive5000.filt.q50.vcf.gz freemuxlet_5000

# Run scSplit
./run_scsplit.sh ../40k_all.bam ../donor_all.bcs varsNaive.filt.q100.vcf.gz scsplit_all
./run_scsplit.sh ../40k_sc_2000.chosen.bam ../samp2000.barcodes2 varsNaive2000.filt.q50.vcf.gz scsplit_2000
./run_scsplit.sh ../40k_sc_5000.chosen.bam ../samp5000.barcodes2 varsNaive5000.filt.q50.vcf.gz scsplit_5000

# Run souporcell
./run_souporcell.sh ../40k_all.bam ../donor_all.bcs souporcell_all
./run_souporcell.sh ../40k_sc_2000.chosen.bam ../samp2000.barcodes2 souporcell_2000
./run_souporcell.sh ../40k_sc_5000.chosen.bam ../samp5000.barcodes2 souporcell_5000


# Run Vireo
./run_vireo.sh ../40k_all.bam ../donor_all.bcs varsNaive.filt.q100.vcf.gz vireo_all
./run_vireo.sh ../40k_sc_2000.chosen.bam ../samp2000.barcodes2 varsNaive2000.filt.q50.vcf.gz vireo_2000
./run_vireo.sh ../40k_sc_5000.chosen.bam ../samp5000.barcodes2 varsNaive5000.filt.q50.vcf.gz vireo_5000

# Get correct answers
cat ../ids_true_consensus.txt | grep -v -F "+" | awk '{printf("%s\t%s\tS\t100\n", $1, $2);}' > ids_true.assignments
cat ../ids_true_consensus.txt | grep -F "+" | awk '{printf("%s\t%s\tD\t100\n", $1, $2);}' >> ids_true.assignments

grep -w -f ../samp2000.barcodes ids_true.assignments > ids_true_2000.assignments
grep -w -f ../samp5000.barcodes ids_true.assignments > ids_true_5000.assignments

# Run all comparisons
${CELLBOUNCER}/plot/compare_assignments.R ids_true.assignments MTvars_all demux_vcf_all_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_2000.assignments mitovars_2000 demux_vcf_2000_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_5000.assignments mitovars_5000 demux_vcf_5000_comp D

${CELLBOUNCER}/plot/compare_assignments.R ids_true.assignments MTvars_all_refined demux_vcf_all_refined_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_2000.assignments mitovars_2000_refined demux_vcf_2000_refined_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_5000.assignments mitovars_5000_refined demux_vcf_5000_refined_comp D

${CELLBOUNCER}/plot/compare_assignments.R ids_true.assignments freemuxlet_all freemuxlet_all_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_2000.assignments freemuxlet_2000 freemuxlet_2000_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_5000.assignments freemuxlet_5000 freemuxlet_5000_comp D

${CELLBOUNCER}/plot/compare_assignments.R ids_true.assignments scsplit_all scsplit_all_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_2000.assignments scsplit_2000 scsplit_2000_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_5000.assignments scsplit_5000 scsplit_5000_comp D

${CELLBOUNCER}/plot/compare_assignments.R ids_true.assignments souporcell_all souporcell_all_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_2000.assignments souporcell_2000 souporcell_2000_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_5000.assignments souporcell_5000 souporcell_5000_comp D

${CELLBOUNCER}/plot/compare_assignments.R ids_true.assignments vireo_all vireo_all_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_2000.assignments vireo_2000 vireo_2000_comp D
${CELLBOUNCER}/plot/compare_assignments.R ids_true_5000.assignments vireo_5000 vireo_5000_comp D
