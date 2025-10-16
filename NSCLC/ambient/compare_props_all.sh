#! /usr/bin/env bash

echo -e "Source\tFraction\tdonor_1\tdonor_2\tdonor_3\tdonor_4\tdonor_5\tdonor_6\tdonor_7"
fracs=( "0.05" "0.1" "0.15" "0.2" "0.25" )
for d in $(seq 1 7); do
    for f in "${fracs[@]}"; do
        row=$( ${CELLBOUNCER}/utils/compare_props.R ../demux_vcf_sc_1000.bulkprops "d${d}c${f}.bulkprops" | cut -f5 | tail -n +2 | tr '\n' '\t' )
        echo -e "donor_${d}\t${f}\t${row}" 
    done
done
