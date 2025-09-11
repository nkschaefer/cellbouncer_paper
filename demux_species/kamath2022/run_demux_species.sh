#! /usr/bin/env bash

if [ $# -lt 3 ]; then
    >&2 echo "Arguments:"
    >&2 echo "k (k-mer length)"
    >&2 echo "Number of k-mers to sample (-1 for all)"
    >&2 echo "Path to cellranger (unpacked)"
    exit 1
fi

# Build reference
#${CELLBOUNCER}/utils/demux_species_ref.py -f Rattus_norvegicus.GRCr8.cdna.all.fa Macaca_mulatta.Mmul_10.cdna.all.fa Tupaia_belangeri.TREESHREW.cdna.all.fa -n Rat Rhesus Tree_shrew -O -k $1 -N $2 -o RMT${1}_${2}

wl="${3}/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"

# Create YAML file
echo "kmers: \"RMT${1}_${2}\"" > demux_species_${1}_${2}.yml
echo "whitelist: \"${wl}\"" >> demux_species_${1}_${2}.yml
echo "libs: \"libs.txt\"" >> demux_species_${1}_${2}.yml
echo "rna_dir: \"reads\"" >> demux_species_${1}_${2}.yml
echo "output_directory: \"demux_species_${1}_${2}\"" >> demux_species_${1}_${2}.yml
echo "append_libname: false" >> demux_species_${1}_${2}.yml
echo "demux_reads: false" >> demux_species_${1}_${2}.yml
echo "threads: 1" >> demux_species_${1}_${2}.yml

nextflow ${CELLBOUNCER}/pipelines/demux_species.nf -params-file demux_species_${1}_${2}.yml


