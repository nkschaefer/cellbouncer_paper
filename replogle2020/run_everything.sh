#! /usr/bin/env bash

# Data should already be unpacked

# Split relevant experiment up by library
${CELLBOUNCER}/utils/split_mex_libs.py -m GSM4367984_exp6.matrix.mtx.gz -f GSM4367984_exp6.features.tsv.gz -b GSM4367984_exp6.barcodes.tsv.gz -o replogle_split

ds=( "0.01" "0.05" "0.1" "0.25" "0.5" )

# Run on full-count matrices
for idx in $(seq 1 6); do
    # Run CellBouncer (and build .counts file)
    ${CELLBOUNCER}/demux_tags -o replogle_ipsc_split_${idx} --sgRNA -B replogle_split_${idx}.barcodes.tsv.gz \
        -F replogle_split_${idx}.features.tsv.gz -M replogle_split_${idx}.matrix.mtx.gz \
        -t "CRISPR Guide Capture" --libname $idx

    # Run mixture model
    ./pg.py replogle_ipsc_split_${idx}.counts
done

for d in "${ds[@]}"; do
    for idx in $(seq 1 6); do
        ./downsample_mtx.py replogle_split_${idx}.matrix.mtx.gz replogle_split_${idx}.features.tsv.gz $d
        
        # Run CellBouncer (and build .counts file)
        ${CELLBOUNCER}/demux_tags -o replogle_ipsc_split_${idx}_ds${d} --sgRNA -B replogle_split_${idx}.barcodes.tsv.gz \
            -F replogle_split_${idx}.features.tsv.gz -M replogle_split_${idx}_ds${d}.matrix.mtx.gz \
            -t "CRISPR Guide Capture" --libname $idx
        # Run mixture model
        ./pg.py replogle_ipsc_split_${idx}_ds${d}.counts 
    
    done
done

# Merge tables back together
${CELLBOUNCER}/utils/combine_sgrna_tables.py replogle_split_1.table replogle_split_2.table replogle_split_3.table replogle_split_4.table replogle_split_5.table replogle_split_6.table > replogle_combined.table

${CELLBOUNCER}/utils/combine_sgrna_tables.py replogle_split_1.pg.table replogle_split_2.pg.table replogle_split_3.pg.table replogle_split_4.pg.table replogle_split_5.pg.table replogle_split_6.pg.table > replogle_combined.pg.table

for d in "${ds[@]}"; do
    ${CELLBOUNCER}/utils/combine_sgrna_tables.py replogle_split_1_ds${d}.table replogle_split_2_ds${d}.table replogle_split_3_ds${d}.table replogle_split_4_ds${d}.table replogle_split_5_ds${d}.table replogle_split_6_ds${d}.table > replogle_combined_ds${d}.table
    ${CELLBOUNCER}/utils/combine_sgrna_tables.py replogle_split_1_ds${d}.pg.table replogle_split_2_ds${d}.pg.table replogle_split_3_ds${d}.pg.table replogle_split_4_ds${d}.pg.table replogle_split_5_ds${d}.pg.table replogle_split_6_ds${d}.pg.table > replogle_combined_ds${d}.pg.table
done


