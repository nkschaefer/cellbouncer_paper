#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "Please provide the path to a CellRanger human reference (unzipped)"
    exit 1
fi

ref=$1

# Get data
wget https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/067/E-MTAB-9067/Files/E-MTAB-9067.sdrf.txt

# Download data
./download_seqs.sh

nlines=$( cat rnaseq_ids_bc_cellranger.txt | wc -l )

# Add barcodes to reads
seq 1 $nlines | sed 's/^/.\/addbc.sh /' | parallel -j 12 

# Merge samples
mkdir reads_cr
seq 1 15 | while read samp; do
    ./merge_samp.sh $samp
    mv Sample${samp}_bc/Ranzoni*.fastq.gz reads_cr
done

git clone git@github.com:nkschaefer/bio_misc.git
cd bio_misc
make
cd ..

if [ -e bamlist.txt ]; then
    rm bamlist.txt
fi

# Align data with STAR
seq 1 15 | while read samp; do
    ./map.sh ${ref} ${samp}
    ./merge_star.sh $samp
    ./rg.sh $samp
    echo "Sample${samp}_merged_rg.bam" >> bamlist.txt
done

# Call variants - produces VCF
./freebayes.sh $ref

# Create "true" assignments file for VCF refinement
cat cellranger_bc2id.tsv | awk '{printf("%s\t%s\tS\t100", $1, $2);}' > cellranger_bc2id.assignments

