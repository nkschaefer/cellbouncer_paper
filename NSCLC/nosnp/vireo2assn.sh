#! /usr/bin/env bash

cat $1 | tail -n +2 | awk '{if ($2 == "doublet"){ printf("%s\t%s\tD\t%f\n", $1, $7, $8 ); } else if ($2 != "unassigned"){ printf("%s\t%s\tS\t%f\n", $1, $2, $3); } }' | sed 's/-1//' | sed 's/,/\t/' | awk '{if (NF == 5){ if ($2 < $3){ printf("%s\t%s+%s\t%s\t%s\n", $1, $2, $3, $4, $5); } else{ printf("%s\t%s+%s\t%s\t%s\n", $1, $3, $2, $4, $5); } } else{ print $0; }}'

