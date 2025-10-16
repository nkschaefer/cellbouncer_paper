#! /usr/bin/env bash

cat $1 | cut -f2,5,6,10 | tail -n +2 | sed 's/,/\t/g' | sed 's/-1//' | awk '{if ($2 == "DBL"){ if ($3 < $4){ printf("%s\t%s+%s\tD\t%f\n", $1, $3, $4, $6); } else{ printf("%s\t%s+%s\tD\t%f\n", $1, $4, $3, $6); }} else{ printf("%s\t%s\tS\t%f\n", $1, $3, $6); }}'
