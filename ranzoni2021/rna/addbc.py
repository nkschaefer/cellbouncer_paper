#! /usr/bin/env python3
import sys
import os
import gzip
import argparse
import random

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--read1", "-1", help="Forward read file", required=True)
    parser.add_argument("--read2", "-2", help="Reverse read file", required=True)
    parser.add_argument("--barcode", "-b", help="Barcode sequence", required=True)
    parser.add_argument("--output_prefix", "-o", help="Output base file name", required=True)
    return parser.parse_args()

def next_read(fn):
    f = gzip.open(fn, 'r')
    read = {'id': None, 'seq': None, 'qual': None}
    for line in f:
        line = line.decode().rstrip()
        if line != '+':
            if read['id'] is None:
                read['id'] = line.split(' ')[0].split('@')[1]
            elif read['seq'] is None:
                read['seq'] = line
            elif read['qual'] is None:
                read['qual'] = line
                yield read
                read['id'] = None
                read['seq'] = None
                read['qual'] = None

def main(args):
    options = parse_args()
    
    out1 = gzip.open('{}_R1.fastq.gz'.format(options.output_prefix), 'wt')
    out2 = gzip.open('{}_R2.fastq.gz'.format(options.output_prefix), 'wt')

    gen1 = next_read(options.read1)
    gen2 = next_read(options.read2)

    bases = ['A', 'C', 'G', 'T']
    R1qual = "?" * 50

    num = 0
    for seq1 in gen1:
        seq2 = next(gen2)
        
        R1seq1 = list(options.barcode)
        for i in range(0,34):
            R1seq1.append(random.choice(bases))
        R1seq2 = list(options.barcode)
        for i in range(0,34):
            R1seq2.append(random.choice(bases))
        
        R1seq1 = "".join(R1seq1)
        R1seq2 = "".join(R1seq2)

        # Create barcode read
        print("@{}.1".format(seq1['id']), file=out1)
        print(R1seq1, file=out1)
        print("+", file=out1)
        print(R1qual, file=out1)

        print("@{}.2".format(seq1['id']), file=out1)
        print(R1seq2, file=out1)
        print("+", file=out1)
        print(R1qual, file=out1)
        
        print("@{}.1".format(seq1['id']), file=out2)
        print(seq1['seq'], file=out2)
        print("+", file=out2)
        print(seq1['qual'], file=out2)

        print("@{}.2".format(seq2['id']), file=out2)
        print(seq2['seq'], file=out2)
        print("+", file=out2)
        print(seq2['qual'], file=out2)

    out1.close()
    out2.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
