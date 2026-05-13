#! /usr/bin/env bash

cat combined_20_10000000/key.txt | tail -n +2 | grep -v -F "+" | awk '{ printf("%s\t%s\tS\t100\n", $1, $2);}' > true_20_10m.assignments
cat combined_20_10000000/key.txt | tail -n +2 | grep -F "+" | awk '{ printf("%s\t%s\tD\t100\n", $1, $2);}' >> true_20_10m.assignments

${CELLBOUNCER}/plot/compare_assignments.R true_20_10m.assignments combined_20_10000000/species.filt.assignments k20_10m_comp D

