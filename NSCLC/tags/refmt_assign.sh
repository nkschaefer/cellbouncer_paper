#! /usr/bin/env bash
cat $1 | awk '{ if ($3 != "S"){ printf("%s\tmultiplet\n", $1); } else{ printf("%s\t%s\n", $1, $2); }}'

