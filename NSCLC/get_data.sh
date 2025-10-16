#! /usr/bin/env bash

# Create one merged BAM
bams=$( ls "40k_NSCLC_DTC_3p_HT_nextgem_donor_*_count_sample_alignments.bam" )
samtools merge 40k_all.bam $bams
samtools index 40k_all.bam

# Add sample read groups to BAM files
git clone git@github.com:nkschaefer/bio_misc.git
cd bio_misc
make
cd ..
if [ -e "bams.txt" ]; then
    rm bams.txt
fi
for donor in $(seq 1 7); do
    bio_misc/bam_dummy_rg -b 40k_NSCLC_DTC_3p_HT_nextgem_donor_${donor}_count_sample_alignments.bam -s donor_${donor} -l donor_${donor} -o donor_${donor}.bam
    samtools index donor_${donor}.bam
    echo "donor_${donor}.bam" >> bams.txt
done

# Get reference genome
wget https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
samtools faidx hg38.fa

# Call variants with freebayes (12 parallel processes)
cat hg38.fa.fai | cut -f1-2 | grep -v chrUn | grep -v random | grep -v chrM | grep -v chrY | grep -v alt | while read line; do
    chrom=$( echo -e "$line" | cut -f1 )
    chromlen=$( echo -e "$line" | cut -f2 )
    echo "freebayes -f hg38.fa -L bams.txt -r ${chrom}:0-${chromlen} > vars.${chrom}.vcf" 
done | parallel -j 12

# Compile variants
cat vars.chr21.vcf | grep -P "^#" > vars.concat.vcf
seq 1 22 | while read chrom; do
    cat vars.chr${chrom}.vcf | grep -v -P "^#" >> vars.concat.vcf
done
cat vars.chrX.vcf | grep -v -P "^#" >> vars.concat.vcf
bgzip vars.concat.vcf
tabix -p vcf vars.concat.vcf.gz

# Filter variants
bcftools view -c 1 -C 13 -v snps -i 'QUAL > 100 & F_MISSING < 0.5' -O z vars.concat.vcf.gz > 500k.vcf.gz
tabix -p vcf 500k.vcf.gz

# Downsample variants
cat vars.chr21.vcf | grep -P "^#" > 250k.vcf
cp 250k.vcf 100k.vcf
cp 250k.vcf 50k.vcf
cp 250k.vcf 25k.vcf
cp 250k.vcf 10k.vcf
cp 250k.vcf 5k.vcf
cp 250k.vcf 1k.vcf
zcat 500k.vcf.gz | grep -v -P "^#" | sort -R | head -n 250000 | sort -k1,1V -k2,2n >> 250k.vcf
zcat 500k.vcf.gz | grep -v -P "^#" | sort -R | head -n 100000 | sort -k1,1V -k2,2n >> 100k.vcf
zcat 500k.vcf.gz | grep -v -P "^#" | sort -R | head -n 50000 | sort -k1,1V -k2,2n >> 50k.vcf
zcat 500k.vcf.gz | grep -v -P "^#" | sort -R | head -n 25000 | sort -k1,1V -k2,2n >> 25k.vcf
zcat 500k.vcf.gz | grep -v -P "^#" | sort -R | head -n 10000 | sort -k1,1V -k2,2n >> 10k.vcf
zcat 500k.vcf.gz | grep -v -P "^#" | sort -R | head -n 5000 | sort -k1,1V -k2,2n >> 5k.vcf
zcat 500k.vcf.gz | grep -v -P "^#" | sort -R | head -n 1000 | sort -k1,1V -k2,2n >> 1k.vcf
bgzip 250k.vcf
bgzip 100k.vcf
bgzip 50k.vcf
bgzip 25k.vcf
bgzip 10k.vcf
bgzip 5k.vcf
bgzip 1k.vcf
tabix -p vcf 250k.vcf.gz
tabix -p vcf 100k.vcf.gz
tabix -p vcf 50k.vcf.gz
tabix -p vcf 25k.vcf.gz
tabix -p vcf 10k.vcf.gz
tabix -p vcf 5k.vcf.gz
tabix -p vcf 1k.vcf.gz

# Create downsampled BAMs
samtools view -s 0.75 -bh 40k_all.bam 40k_s0.75.bam
samtools view -s 0.5 -bh 40k_all.bam 40k_s0.5.bam
samtools view -s 0.25 -bh 40k_all.bam 40k_s0.25.bam
samtools view -s 0.1 -bh 40k_all.bam 40k_s0.1.bam
samtools index 40k_s0.75.bam
samtools index 40k_s0.5.bam
samtools index 40k_s0.25.bam
samtools index 40k_s0.1.bam

# Create downsampled sets of cells
cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_assignment_confidence_table.csv | sed 's/,/\t/g' | cut -f9,12 | tail -n +2 | grep -v Blanks | grep -v Unassigned | sed 's/-1//' > barcodes.tsv
cat barcodes.tsv | sort -R | head -100 | awk '{printf("%s\tchosen\tS\t100\n");}' > samp100.assignments
cat barcodes.tsv | sort -R | head -500 | awk '{printf("%s\tchosen\tS\t500\n");}' > samp500.assignments
cat barcodes.tsv | sort -R | head -1000 | awk '{printf("%s\tchosen\tS\t1000\n");}' > samp1000.assignments
cat barcodes.tsv | sort -R | head -2000 | awk '{printf("%s\tchosen\tS\t2000\n");}' > samp2000.assignments
cat barcodes.tsv | sort -R | head -5000 | awk '{printf("%s\tchosen\tS\t5000\n");}' > samp5000.assignments
cat barcodes.tsv | sort -R | head -10000 | awk '{printf("%s\tchosen\tS\t10000\n");}' > samp10000.assignments
cat barcodes.tsv | sort -R | head -20000 | awk '{printf("%s\tchosen\tS\t20000\n");}' > samp20000.assignments
cat samp100.assignments | cut -f1 > samp100.barcodes
cat samp500.assignments | cut -f1 > samp500.barcodes
cat samp1000.assignments | cut -f1 > samp1000.barcodes
cat samp2000.assignments | cut -f1 > samp2000.barcodes
cat samp5000.assignments | cut -f1 > samp5000.barcodes
cat samp10000.assignments | cut -f1 > samp10000.barcodes
cat samp20000.assignments | cut -f1 > samp20000.barcodes

# Make versions of barcodes files that have the -1 appended to each
cat samp100.barcodes | sed 's/$/-1/' > samp100.barcodes2
cat samp500.barcodes | sed 's/$/-1/' > samp500.barcodes2
cat samp1000.barcodes | sed 's/$/-1/' > samp1000.barcodes2
cat samp2000.barcodes | sed 's/$/-1/' > samp2000.barcodes2
cat samp5000.barcodes | sed 's/$/-1/' > samp5000.barcodes2
cat samp10000.barcodes | sed 's/$/-1/' > samp10000.barcodes2
cat samp20000.barcodes | sed 's/$/-1/' > samp20000.barcodes2

${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp100.assignments -o 40k_sc_100
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp500.assignments -o 40k_sc_500
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp1000.assignments -o 40k_sc_1000
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp2000.assignments -o 40k_sc_2000
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp5000.assignments -o 40k_sc_5000
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp10000.assignments -o 40k_sc_10000
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp20000.assignments -o 40k_sc_20000
samtools index 40k_sc_100.chosen.bam
samtools index 40k_sc_500.chosen.bam
samtools index 40k_sc_1000.chosen.bam
samtools index 40k_sc_2000.chosen.bam
samtools index 40k_sc_5000.chosen.bam
samtools index 40k_sc_10000.chosen.bam
samtools index 40k_sc_20000.chosen.bam

cat 40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_multiplexing_analysis_tag_calls_per_cell.csv | tail -n +2 | cut -d',' -f1 > donor_all.bcs

# Get mitochondria-specific BAM and extract "junk" droplet barcodes from it
samtools view -bh 40k_all.bam chrM > chrM.bam
samtools index chrM.bam
samtools view chrM.bam | grep -P -oh "CB:Z:[ACGT]+-1" | cut -d':' -f3 | sort | uniq | grep -v -f donor_all.bcs > junk.bcs

# Create "plus junk" barcode lists for mitochondrial demultiplexing
cat samp2000.assignments > samp2000.plusjunk.assignments
cat junk.bcs | sort -R | head -49890 | awk '{ printf("%s\tchosen\tS\t100\n", $1); }' >> samp2000.plusjunk.assignments
cat samp5000.assignments > samp5000.plusjunk.assignments
cat junk.bcs | sort -R | head -114934 | awk '{ printf("%s\tchosen\tS\t100\n", $1); }' >> samp5000.plusjunk.assignments

# Create BAM files
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp2000.plusjunk.assignments -o 40k_sc_2000.plusjunk
${CELLBOUNCER}/utils/bam_split_bcs -b 40k_all.bam -a samp5000.plusjunk.assignments -o 40k_sc_5000.plusjunk
samtools index 40k_sc_2000.plusjunk.chosen.bam
samtools index 40k_sc_5000.plusjunk.chosen.bam

