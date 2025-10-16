#! /usr/bin/env bash

cat E-MTAB-9068.sdrf.txt | tail -n +2 | while read line; do
    id=$( echo -e "$line" | cut -f2 )
    sample=$( echo -e "$line" | cut -f12 )
    url=$( echo -e "$line" | cut -f39 )
    if [ ! -d ${sample} ]; then
        mkdir ${sample}
    fi
    cd ${sample}
    wget $url
    cd ..
done
