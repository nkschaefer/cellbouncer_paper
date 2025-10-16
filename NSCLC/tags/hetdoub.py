#! /usr/bin/env python3
import sys
import os

def main(args):
    tot = 0
    m = 0
    for line in sys.stdin:
        line = line.rstrip()
        if line != "":
            bc, assn = line.split('\t')
            if assn != "negative":
                tot += 1
                if assn == "multiplet":
                    m += 1
    print(m/tot)

if __name__ == '__main__':
    sys.exit(main(sys.argv))
