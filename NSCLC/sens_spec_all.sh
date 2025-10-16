#! /usr/bin/env bash

if [ -e "sens_spec_ncell.txt" ]; then
    rm sens_spec_ncell.txt
fi

if  [ -e "sens_spec_downsample.txt" ]; then
    rm sens_spec_downsample.txt
fi

if [ -e "sens_spec_nvar.txt" ]; then
    rm sens_spec_nvar.txt
fi

# ncells
ncell=( 100 500 1000 5000 10000 20000 )
for n in ${ncell[@]}; do
    >&2 echo "ncell $n"
    ./sens_spec.R demuxlet_sc_${n}.assignments >> sens_spec_ncell.txt
    ./sens_spec.R demuxalot_sc_${n}.assignments >> sens_spec_ncell.txt
    ./sens_spec.R vireo_sc_${n}.assignments >> sens_spec_ncell.txt
    ./sens_spec.R demux_vcf_sc_${n}.assignments >> sens_spec_ncell.txt
    ./sens_spec.R demux_vcf2_sc_${n}.decontam.assignments >> sens_spec_ncell.txt
done

# nvars
nvar=( "1k" "10k" "25k" "50k" "100k" "250k" "500k" )
for n in ${nvar[@]}; do
    >&2 echo "nvar $n"
    ./sens_spec.R demuxlet_donor_all_${n}.assignments >> sens_spec_nvar.txt
    ./sens_spec.R demuxalot_donor_all_${n}.assignments >> sens_spec_nvar.txt
    ./sens_spec.R vireo_donor_all_${n}.assignments >> sens_spec_nvar.txt
    ./sens_spec.R demux_vcf_donor_all_${n}.assignments >> sens_spec_nvar.txt
    ./sens_spec.R demux_vcf2_donor_all_${n}.decontam.assignments >> sens_spec_nvar.txt
done

# downsample
ds=( "0.1" "0.25" "0.5" "0.75" )
for n in ${ds[@]}; do
    >&2 echo "ds $n"
    ./sens_spec.R demuxlet_ds${n}.assignments >> sens_spec_downsample.txt
    ./sens_spec.R demuxalot_ds${n}.assignments >> sens_spec_downsample.txt
    ./sens_spec.R vireo_ds${n}.assignments >> sens_spec_downsample.txt
    ./sens_spec.R demux_vcf_ds${n}.assignments >> sens_spec_downsample.txt
    ./sens_spec.R demux_vcf2_ds${n}.decontam.assignments >> sens_spec_downsample.txt
done

