#! /usr/bin/env bash

if [ $# -lt 3 ]; then
    >&2 echo "Required arguments:"
    >&2 echo "path to (unpacked) human CellRanger ref data"
    >&2 echo "path to (unpacked) mouse CellRanger ref data"
    >&2 echo "path to (unpacked) CellRanger programs"
    exit 1
fi

human=$1
mouse=$2
cellranger=$3

# Build composite reference genome (FASTA)
cat ${human}/fasta/genome.fa | sed 's/>/>hg38_/' > mouse_human_2020.fa
cat ${mouse}/fasta/genome.fa | sed 's/>/>mm10_/' >> mouse_human_2020.fa

# Build composite annotation (GTF)
cat ${human}/genes/genes.gtf | grep -v "#" | sed 's/^/hg38_/' > mouse_human_2020.gtf
cat ${mouse}/genes/genes.gtf | grep -v "#" | sed 's/^/mm10_/' >> mouse_human_2020.gtf

# Get barcode whitelist
zcat ${cellranger}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz > wl.txt

# Build STAR index
./star_index.sh

# Map data with STARsolo
./starsolo.sh

# Extract data to run demux_species
bam="mouse_human_STARAligned.sortedByCoord.out.bam"

# Run CellBouncer on extracted counts
${CELLBOUNCER}/utils/composite_bam2counts -b $bam -o demux_species_composite
${CELLBOUNCER}/demux_species -d -o demux_species_composite

