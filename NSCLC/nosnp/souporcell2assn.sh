#! /usr/bin/env bash

if [ $# -lt 1 ]; then
    >&2 echo "ERROR: provide souporcell out dir"
    exit 1
fi

cat ${1}/clusters.tsv | sed 's/-1//' | cut -f1-5 | tail -n +2 | grep -v unassigned | awk '{if ($2 == "singlet"){ printf("%s\t%s\tS\t%f\n", $1, $3, $4); } else{ printf("%s\t%s\tD\t%f\n", $1, $3, $5); }}' | sed 's/\//\t/' | awk '{ if ($4 == "D"){ if ($2 < $3){ printf("%s\t%s+%s\t%s\t%f\n", $1, $2, $3, $4, $5); } else{ printf("%s\t%s+%s\t%s\t%f\n", $1, $3, $2, $4, $5); }} else{ printf("%s\n", $0); }}'
