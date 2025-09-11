#! /usr/bin/env bash

# Consolidate data, make assignments, and hierarchically cluster cells
# for each run

./make_tables.R 15_-1
./make_tables.R 20_-1
./make_tables.R 25_-1
./make_tables.R 35_-1
./make_tables.R 45_-1
./make_tables.R 20_1000000
./make_tables.R 20_5000000
./make_tables.R 20_10000000
./make_tables.R 20_20000000

${CELLBOUNCER}/demux_species -d -o combined_15_-1
${CELLBOUNCER}/plot/species.R combined_15_-1
${CELLBOUNCER}/demux_species -d -o combined_20_-1
${CELLBOUNCER}/plot/species.R combined_20_-1
${CELLBOUNCER}/demux_species -d -o combined_25_-1
${CELLBOUNCER}/plot/species.R combined_25_-1
${CELLBOUNCER}/demux_species -d -o combined_35_-1
${CELLBOUNCER}/plot/species.R combined_35_-1
${CELLBOUNCER}/demux_species -d -o combined_45_-1
${CELLBOUNCER}/plot/species.R combined_45_-1
${CELLBOUNCER}/demux_species -d -o combined_20_1000000
${CELLBOUNCER}/plot/species.R combined_20_1000000
${CELLBOUNCER}/demux_species -d -o combined_20_5000000
${CELLBOUNCER}/plot/species.R combined_20_5000000
${CELLBOUNCER}/demux_species -d -o combined_20_10000000
${CELLBOUNCER}/plot/species.R combined_20_10000000
${CELLBOUNCER}/demux_species -d -o combined_20_20000000
${CELLBOUNCER}/plot/species.R combined_20_20000000

# Now compare these results to truth, and output these data as
# tables
./compare_kmers.R

# Plot overlap between true species IDs and Cellbouncer-filtered IDs
# using 10 million sampled 20-mers
./make_jaccard_plot.R
