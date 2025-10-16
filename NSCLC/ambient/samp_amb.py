#! /usr/bin/env python3
import sys
import os
import subprocess
import gzip
import random
import pysam
import argparse
import glob

def samp_bam(percamb, bam_in, bam_out, totreads, nbam):
    samp_probs = {}
    f = open('samp1000_bccounts.txt', 'r')
    samp_prob_sum = 0.0
    for line in f:
        line = line.rstrip()
        dat = line.strip().split(" ")
        nreads = int(dat[0])
        bc = dat[1]
        
        # How many new reads to sample?
        # x / (n + x) = percamb
        # (n+x)percamb = x
        # np + xp = x
        # np = x - xp
        # np = x(1-p)
        # x = np/(1-p)

        nsamp = (nreads*percamb)/(1-percamb)
        nsamp *= (1/nbam)
        samp_prob = nsamp/totreads
        samp_probs[bc] = samp_prob
        samp_prob_sum += samp_prob

    # Now sort probs 
    intervals = []
    bcs = []

    probsort = []
    for bc in samp_probs:
        probsort.append((samp_probs[bc], bc))
    probsort.sort()
    
    startprob = 0
    for tup in probsort:
        prob, bc = tup
        endprob = startprob + prob
        intervals.append((startprob, endprob))
        startprob = endprob
        bcs.append(bc)
    
    bamf = pysam.AlignmentFile(bam_in, "rb")
    outf = pysam.AlignmentFile(bam_out, "wb", template=bamf)

    for rec in bamf:
        # Skip reads missing barcodes or not meeting criteria
        if rec.is_duplicate or not rec.is_mapped or rec.is_secondary or not rec.has_tag('CB'):
            continue
        r = random.random()
        if r >= startprob:
            # Out of range of thresholds
            continue
        
        ival_found = False
        for i, ival in enumerate(intervals):
            if r >= ival[0] and r < ival[1]:
                # In interval.
                rec.set_tag("CB", bcs[i])
                outf.write(rec)
                ival_found = True
                break

    bamf.close()
    outf.close()

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--amb", "-a", type=float, help="Percent ambient RNA per cell", required=True)
    parser.add_argument("--donors", "-d", nargs="+", help="Donors to include", required=True)
    parser.add_argument("--out", "-o", help="Output prefix", required=True)
    parser.add_argument("--cellbouncer", "-c", help="Path to CellBouncer", required=True)
    return parser.parse_args()

def main(args):
    
    options = parse_args()

    percamb = options.amb
    if percamb <= 0 or percamb >= 1:
        print("ERROR: perc amb must be between 0 and 1", file=sys.stderr)
        exit(1)
    
    # Used samtools view -F1284 -c [BAM]
    # donor 1
    totreads1 = 9224269
    # donor 2
    totreads2 = 19095977
    # donor 3
    totreads3 = 3953362
    # donor 4
    totreads4 = 3950828
    # donor 5
    totreads5 = 13797308
    # donor 7
    totreads7 = 5381775
    # donor 6
    totreads6 = 21913975
    
    nreads = [totreads1, totreads2, totreads3, totreads4, totreads5, totreads6, totreads7]
    bams = []
    bams_out = []
    for donor in range(1, 8):
        bams.append('donor{}_nosamp1000.chosen.vars.bam'.format(donor))
        bams_out.append('donor{}_{}.bam'.format(donor, options.out))
    
    to_merge = []
    ndonor = len(options.donors)
    for donor in options.donors:
        donor_idx = int(donor) - 1
        print("Sample {}...".format(bams[donor_idx]), file=sys.stderr)
        samp_bam(percamb, bams[donor_idx], bams_out[donor_idx], nreads[donor_idx], ndonor)
        to_merge.append(bams_out[donor_idx])
    # Merge 
    print("Merging to {}_merged.bam...".format(options.out), file=sys.stderr)
    subprocess.call(['samtools', 'merge', '-f', '{}_merged.bam'.format(options.out), \
        '40k_sc_1000.chosen.vars.bam'] + to_merge)
    # Index
    subprocess.call(['samtools', 'index', '{}_merged.bam'.format(options.out)])
    
    # Clean up
    for fn in glob.glob('{}.*'.format(options.out)):
        os.unlink(fn)

    #print("Counting variants...", file=sys.stderr)
    #subprocess.call(['{}/demux_vcf'.format(options.cellbouncer), '-v', '../500k.vcf.gz', '-b', \
    #    '{}_merged.bam'.format(options.out), '-o', options.out])
    #subprocess.call(['cp', 'demux_vcf_sc_1000.assignments', '{}.assignments'.format(options.out)])


if __name__ == '__main__':
    sys.exit(main(sys.argv))
