#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: please provide path to CellRanger barcode whitelist (3M-february-2018.txt.gz)"
    >&2 echo "This is in [cellranger]/lib/python/cellranger/barcodes"
    exit 1
fi

./subsamp_all.sh $1

if [ -e downsample.stats ]; then
    rm downsample.stats
fi
if [ -e all_hetdoub.stats ]; then
    rm all_hetdoub.stats
fi 
./run_downsampled.sh full.counts >> downsample.stats
./run_hetdoub.sh full.counts >> all_hetdoub.stats

# Now run all programs on downsampled tag counts
ls s*.counts | while read fn; do
    ./run_downsampled.sh $fn >> downsample.stats
    ./run_hetdoub.sh $fn >> all_hetdoub.stats
done

# Plot downsampled
./plot_downsample.R downsample.stats

if [ -e ambsamp.stats ]; then
    rm ambsamp.stats
fi

if [ -e all_hetdoub_ambsamp.stats ]; then
    rm all_hetdoub_ambsamp.stats
fi

ambs=( "0.1" "0.25" "0.5" "0.75" "0.9" )
for amb in "${ambs[@]}"; do
    amb2=$( echo "$amb" | sed 's/./_/' )
    for rep in $(seq 1 3 ); do
        name="ampsamp_${amb2}"
        if [ $rep -eq 2 ]; then
            name="${name}b"
        elif [ $rep -eq 3 ]; then
            name="${name}c"
        fi
        ./samp_amb.R $amb > "${name}.counts"
        ./run_all_ampbsamp.sh ${name}.counts >> ambsamp.stats 
        ./run_all_hetdoub.sh ${name}.counts >> all_hetdoub_ambsamp.stats
    done
done

# Plot results
./plot_stats_ambsamp.R ambsamp.stats

# Plot heterotypic doublet rates
./plot_hetdoub.R 
