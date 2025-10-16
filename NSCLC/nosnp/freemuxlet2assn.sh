#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ARGS: output prefix"
    exit 1
fi

name="${1}_freemuxlet.clust1.samples.gz"
zcat $name | cut -f2,5,6,10 | tail -n +2 | sed 's/,/\t/g' | sed 's/-1//' | awk '{if ($2 == "DBL"){ if ($3 < $4){ printf("%s\t%s+%s\tD\t%f\n", $1, $3, $4, $6); } else{ printf("%s\t%s+%s\tD\t%f\n", $1, $4, $3, $6); }} else{ printf("%s\t%s\tS\t%f\n", $1, $3, $6); }}'
