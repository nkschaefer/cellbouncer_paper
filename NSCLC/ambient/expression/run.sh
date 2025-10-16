#! /usr/bin/env bash
if [ $# -lt 3 ]; then
    >&2 echo "REQUIRED: filtered_feature_bc_matrix directory (unpacked), raw_feature_bc_matrix directory (unpacked) cellranger human refernece (unpacked)"
    exit 1
fi

filt=$1
raw=$2
ref=$3

# Subset raw matrix to make cellbender happy
${CELLBOUNCER}/utils/subs_mex_featuretype.py -m ${2}/matrix.mtx.gz -f ${2}/features.tsv.gz -b ${2}/barcodes.tsv.gz -t "Gene Expression" -o raw_gex

./cellbender.sh raw_gex

# This will produce clusters needed by SoupX and quant_contam
./decontX.R $filt
./soupX.R $filt $raw

# Run quant_contam to remove amb RNA
${CELLBOUNCER}/quant_contam -o ../../demux_vcf_donor_all_500k -T 12 -M ${filt}/matrix.mtx.gz -F ${filt}/features.tsv.gz -B ${filt}/barcodes.tsv.gz -t "Gene Expression" -c decontX_clusters.txt

# Compile data
./compile_expr.R

# Make plot of correlation between contamination rates & gene expr profiles
./contam_corr.R

# Run velocyto
zcat ${ref}/genes/genes.gtf.gz > genes.gtf

### IMPORTANT: edit this to your conda installation of velocyto
source ~/miniforge3/bin/activate velocyto
velocyto run -o velocyto ../../40k_all.bam genes.gtf

# Extract spliced/unspliced per gene from velocyto output
loomfile=$( ls velocyto/*.loom | head -1 )
./extract_velocyto.R $loomfile > gene_splice.txt

# Compile data for lncAtlas plot
./compile_lncatlas.R
# Make plot
./plot_lncatlas.R

# Run Wilcoxon tests of ambient RNA expression vector loadings for genes
# in and not in a previously published list of neuronal ambient RNA genes
./caglayan2022_wilcoxon.R

# Look at "reverse" bulkprops - log likelihood of genotypic makeup at individual
# genes given the consensus bulk proportions genome-wide
${CELLBOUNCER}/bulkprops -p ../../demux_vcf_donor_all_500k.bulkprops -b ../../40k_all.bam -v ../../500k.vcf.gz -e 0.001 -g > bulkprops_ll.txt

# Make plot
./make_ll_plot.R



