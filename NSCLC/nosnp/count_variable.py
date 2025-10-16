#! /usr/bin/env python3
import sys
import os
from collections import Counter

def main(args):
    nvar = 0
    for line in sys.stdin:
        line = line.rstrip()
        if line[0] != "#":
            dat = line.split('\t')
            counts = Counter()
            for fld in dat[9:]:
                gt = fld.split(':')[0]
                counts[gt] += 1
            if len(counts.keys()) > 1:
                nvar += 1
            
    print(nvar)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
