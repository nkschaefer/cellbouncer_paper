#! /usr/bin/env bash

if [ $# -lt 3 ]; then
    >&2 echo "ARGS: R1 subsamp rep"
    exit 1
fi
r1=$1
samp=$2
rep=$3

r2=$( echo -e "$r1" | sed 's/_R1_/_R2_/' )

dirn="${r1%/*}"

# Keep this consistent so reads stay paired
seed=100
dir_out="s${samp}_${dirn}"
if [ "$rep" == "2" ]; then
    dir_out="s${samp}b_${dirn}"
    seed=200
elif [ "$rep" == "3" ]; then
    dir_out="s${samp}c_${dirn}"
    seed=300
fi
r1base="${r1##*/}"
r2base="${r2##*/}"

if [ ! -d $dir_out ]; then
    mkdir $dir_out
fi

source ~/miniforge3/bin/activate seqkit
seqkit sample -j 1 -p $samp -s $seed $r1 -o "${dir_out}/${r1base}"
seqkit sample -j 1 -p $samp -s $seed $r2 -o "${dir_out}/${r2base}"

