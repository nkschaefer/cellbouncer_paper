#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "Please provide the path to CellRanger (unpacked)."
    exit 1
fi

cellranger=$1

# Create yml for running demux_species
echo "kmers: \"HM32_20M\"" > demux_species.yml
echo "whitelist: \"${cellranger}/lib/python/cellranger/barcodes/737K-august-2016.txt\"" >> demux_species.yml
echo "libs: \"libs.txt\"" >> demux_species.yml
echo "rna_dir: \"geo_submission_hellmann_18082023\"" >> demux_species.yml
echo "output_directory: \"demux_species_out\"" >> demux_species.yml
echo "num_chunks: 100" >> demux_species.yml
echo "demux_reads: false" >> demux_species.yml

# Run demux_species
nextflow ${CELLBOUNCER}/pipelines/demux_species.nf -params-file demux_species.yml

# Concatenate data from the two libraries
if [ ! -d demux_species_out/combined ]; then
    mkdir demux_species_out/combined
fi
./combine_runs.R demux_species_out > demux_species_out/combined/species_counts.txt
cd demux_species_out
cp lane1_SI_TT_G1/species_names.txt combined 
${CELLBOUNCER}/demux_species -o combined -d
cd ..

# Process tag counts from Edenhofer 2024
# Reformat the provided table, collapse sgRNAs by whether they target human, macaque or both
# Print counts for cells with >= 10 sgRNA reads to a table in the format recognized by
# demux_species
mkdir edenhofer2024_demux_species
./refmt_sgrna_data.py > edenhofer2024_demux_species
echo -e "0\tHuman" > edenhofer2024_demux_species
echo -e "1\tMacaque" >> edenhofer2024_demux_species
echo -e "2\tBoth" >> edenhofer2024_demux_species

# Run demux_species to assign species identity to Edenhofer 2024 data from tag counts
${CELLBOUNCER}/demux_species -o edenhofer2024_demux_species -d

# Compare results (and plot Jaccard indices)
${CELLBOUNCER}/plot/compare_assignments.R demux_species_out/combined/species.assignments edenhofer2024_demux_species/species.assignments compare_edenhofer2024 D

# Compute precision & recall
./precision_recall.R > precision_recall.txt


