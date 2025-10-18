#! /usr/bin/env python3
import sys
import os
from collections import Counter, defaultdict
import numpy as np

def main(args):
    f = open('flank_corr_dat.tsv', 'r')
    outf = open('flank_corr.tsv', 'w')

    hdr = None
    first = True

    order = ['MT-ND1', 'MT-ND2', 'MT-CO1', 'MT-CO2', 'MT-ATP8', 'MT-ATP6', 'MT-CO3', 'MT-ND3', 'MT-ND4L', 'MT-ND4', 'MT-ND5', 'MT-CYB']
    coords = [(3306, 4262), (4469,5511), (5903,7445), (7585,8269),(8365,8572),(8526,9207),(9206,9990),(10058,10404),\
            (10469,10766), (10759,12137), (12336,14148), (14746,15887)]

    gene2idx = {}
    gene2coords = {}
    for i in range(0, len(order)):
        gene2idx[order[i]] = i
        gene2coords[order[i]] = coords[i]
    
    splits = []
    for i in range(0, len(order)-1):
        j = i + 1
        g1 = order[i]
        g2 = order[j]
        splits.append((g1, g2))

    splitcorrs = {}
    splitcorrs_weights = {}

    for line in f:
        line = line.rstrip()
        dat = line.split('\t')
        if first:
            hdr = {}
            for i in range(0, len(dat)):
                hdr[dat[i]] = i
            first = False
        else:
            for split in splits:
                g1, g2 = split
                g1i = gene2idx[g1]
                g2i = gene2idx[g2]
                
                if dat[hdr['gene1']] not in gene2idx or dat[hdr['gene2']] not in gene2idx:
                    continue

                gene1i = gene2idx[dat[hdr['gene1']]]
                gene2i = gene2idx[dat[hdr['gene2']]]
                
                if (gene1i == g1i and gene2i == g2i) or (gene2i == g1i and gene1i == g2i):
                # It works if the split is within the range
                    for key in hdr:
                        if key != 'gene1' and key != 'gene2' and key != 'dist' and key != 'start' and key != 'end':
                            if key not in splitcorrs:
                                splitcorrs[key] = defaultdict(list)
                                splitcorrs_weights[key] = defaultdict(list)
                            idx = hdr[key]
                            splitcorrs[key][split].append(float(dat[idx]))
    f.close()

    for split in splits:
        for species in splitcorrs:
            meanval = np.mean(np.array(splitcorrs[species][split]))
            w1 = 1
            w2 = 1
            midpt = (gene2coords[split[0]][1] + gene2coords[split[1]][0])/2.0
            print("{}\t{}\t{}\t{}\t{}\t-1\t{}".format(split[0], split[1], midpt, meanval, w1/w2, species), 
                    file=outf)
    
    outf.close()

if __name__ == '__main__':
    sys.exit(main(sys.argv))
