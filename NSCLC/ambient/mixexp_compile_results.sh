#! /usr/bin/env bash

progs=( "demuxlet" "demuxalot" "vireo" )
amounts=( "0.05" "0.1" "0.15" "0.2" "0.25" )
if [ -e "mixexp_compiled.txt" ]; then
    rm mixexp_compiled.txt
fi
for d in $(seq 1 7); do
    for a in ${amounts[@]}; do
        >&2 echo "d${d}c${a}"
        for p in ${progs[@]}; do
            cat "d${d}c${a}_${p}.assignments" | ./sens_spec_1000.R | sed "s/^/${p}\t${a}\tdonor${d}\t/" >> mixexp_compiled.txt
        done
        cat "d${d}c${a}.assignments" | ./sens_spec_1000.R | sed "s/^/CellBouncer\t${a}\tdonor${d}\t/" >> mixexp_compiled.txt
        cat "d${d}c${a}.decontam.assignments" | ./sens_spec_1000.R | sed "s/^/CellBouncer_contam\t${a}\tdonor${d}\t/" >> mixexp_compiled.txt
    done
done
