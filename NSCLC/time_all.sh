#! /usr/bin/env bash

if [ -e "time_ncell.txt" ]; then
    rm time_ncell.txt
fi

if  [ -e "time_downsample.txt" ]; then
    rm time_downsample.txt
fi

if [ -e "time_nvar.txt" ]; then
    rm time_nvar.txt
fi

# ncells
ncell=( 100 500 1000 5000 10000 20000 )
for n in ${ncell[@]}; do
    >&2 echo "ncell $n"
    ./parse_time.py demuxlet_sc_${n}.time >> time_ncell.txt
    ./parse_time.py demuxalot_sc_${n}.time >> time_ncell.txt
    ./parse_time.py vireo_sc_${n}.time >> time_ncell.txt
    ./parse_time.py demux_vcf_sc_${n}.time >> time_ncell.txt
done

# nvars
nvar=( "1k" "10k" "25k" "50k" "100k" "250k" "500k" )
for n in ${nvar[@]}; do
    >&2 echo "nvar $n"
    ./parse_time.py demuxlet_donor_all_${n}.time >> time_nvar.txt
    ./parse_time.py demuxalot_donor_all_${n}.time >> time_nvar.txt
    ./parse_time.py vireo_donor_all_${n}.time >> time_nvar.txt
    ./parse_time.py demux_vcf_donor_all_${n}.time >> time_nvar.txt
done

# downsample
ds=( "0.1" "0.25" "0.5" "0.75" )
for n in ${ds[@]}; do
    >&2 echo "ds $n"
    ./parse_time.py demuxlet_ds${n}.time >> time_downsample.txt
    ./parse_time.py demuxalot_ds${n}.time >> time_downsample.txt
    ./parse_time.py vireo_ds${n}.time >> time_downsample.txt
    ./parse_time.py demux_vcf_ds${n}.time >> time_downsample.txt
done

