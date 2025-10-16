#! /usr/bin/env bash

# Run quant_contam and bulkprops on downsampled set of 1000 cells
${CELLBOUNCER}/quant_contam -o ../demux_vcf_sc_1000
${CELLBOUNCER}/bulkprops -b ../40k_sc_1000.chosen.bam -v ../500k.vcf.gz -e 0.001 -o ../demux_vcf_sc_1000

# Use the downsampled set of 1000 cells

# Get a list of counts of each barcode in the sampled set
cat ../samp1000.barcodes2 | sort | uniq -c > samp1000_bccounts.txt

# For each donor, remove the sampled cells and only keep reads overlapping variant sites.
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | tail -n +2 | grep -w CMO301 | cut -d',' -f9 > donor1.bcs
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | tail -n +2 | grep -w CMO302 | cut -d',' -f9 > donor2.bcs
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | tail -n +2 | grep -w CMO303 | cut -d',' -f9 > donor3.bcs
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | tail -n +2 | grep -w CMO304 | cut -d',' -f9 > donor4.bcs
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | tail -n +2 | grep -w CMO306 | cut -d',' -f9 > donor5.bcs
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | tail -n +2 | grep -w CMO307 | cut -d',' -f9 > donor6.bcs
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | tail -n +2 | grep -w CMO308 | cut -d',' -f9 > donor7.bcs

# Get a BED file of variant sites
zcat ../500k.vcf.gz | grep -v "#" | awk '{printf("%s\t%d\t%d\n", $1, $2-1, $2);}' | bedtools sort -i stdin > 500k.bed

# Subset the chosen 1000 cells to reads hitting variant positions
samtools view -bh --regions-file 500k.bed ../40k_sc_1000.chosen.bam > 40k_sc_1000.chosen.vars.bam
samtools index 40k_sc_1000.chosen.vars.bam

for idx in $(seq 1 7); do
    cat donor${idx}.bcs | grep -v -f ../samp1000.barcodes2 | awk '{printf("%s\tchosen\tS\t100\n", $1); }' > donor${idx}_nosamp1000.assignments
    # Subset the BAM to remove chosen cells
    ${CELLBOUNCER}/utils/bam_split_bcs -b ../40k_NSCLC_DTC_3p_HT_nextgem_donor_${idx}_count_sample_alignments.bam -a donor${idx}_nosamp1000.assignments -o donor${idx}_nosamp1000
    samtools index donor${idx}_nosamp1000.chosen.bam
    
    # Subset the BAM to keep reads hitting variant sites
    samtools view -bh --regions-file 500k.bed donor${idx}_nosamp1000.chosen.bam > donor${idx}_nosamp1000.chosen.vars.bam
    samtools index donor${idx}_nosamp1000.chosen.vars.bam 
done

# Create BAM files with chosen 1000 cells plus simulated ambient RNA
ambs=( "0.05" "0.1" "0.15" "0.2" "0.25" )
for idx in $(seq 1 7 ); do
    for a in "${ambs[@]}"; do
        ./samp_amb.py -c ${CELLBOUNCER} -d donor${idx} -a $a -o "d${idx}c${a}"
    done
done

# Run CellBouncer on each
./make_mixexp_cellbouncer_jobs.sh | parallel -j 12

# Run demuxalot on each
./make_mixexp_demuxalot_jobs.sh | parallel -j 12

# Run demuxlet on each
./make_mixexp_demuxlet_jobs.sh | parallel -j 12

# Run vireo on each
./make_mixexp_vireo_jobs.sh | parallel -j 12

# Run CellBouncer quant_contam on each
./make_quant_contam_jobs.sh | parallel -j 12

# Copy over the "right answer" assignments & run CellBouncer quant_contam
./mixexp_make_corrass.sh

# Run bulkprops on each
./mixexp_bulkprops.sh

# Compile results
./mixexp_compile_results.sh

# Make plot of accuracy of all programs after incorporating ambient RNA knowledge
./plot_mixexp.R

# Make plot of proportions of cells with each assignment & ambient RNA profiles in each
# simulated contaminated data set
# This also creates the plot of inferred contamination rates
./plot_mixexp_props.R

# Generate Table S4 (bulk props significance testing)
./compare_props_all.sh > bulkprops_signif.tsv


