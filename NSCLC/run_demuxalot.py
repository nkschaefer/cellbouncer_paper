#! /usr/bin/env python3
import sys
import os
from demuxalot import Demultiplexer, BarcodeHandler, ProbabilisticGenotypes, count_snps
import gzip
import contextlib

def refmt_df(df):
    """
    Take output from demuxalot and reformat to same format
    as CellBouncer programs, write to stdout.
    """
    df['id'] = df.idxmax(axis=1)
    df['max'] = df.drop(['id'], axis=1).max(axis=1)
    df['second'] = df.drop(['id', 'max'], axis=1).apply(lambda row: row.nlargest(2).values[-1],axis=1)
    df = df.loc[df['max'] != df['second'],:]
    df['type'] = 'S'
    df.loc[df.id.str.contains("+", regex=False),'type'] = 'D'

    df_out = df[['id', 'type', 'max']]
    df_out.to_csv(sys.stdout, sep='\t', header=False, index=True)

def main(args):
    
    if len(args) < 4:
        print("USAGE: run_demuxalot.py bam barcodes VCF", file=sys.stderr)
        exit(1)
    
    bam = args[1]
    bc = args[2]
    vcf = args[3]
    
    samps = ['donor_1', 'donor_2', 'donor_3', 'donor_4', 'donor_5', 'donor_6', 'donor_7']
    barcode_handler = BarcodeHandler.from_file(bc)
    
    genotypes = None
    snps = None
    likelihoods = None
    posterior_probabilities = None
    
    # Redirect garbage to nullfile
    with open(os.devnull, 'w') as devnull:
        with contextlib.redirect_stdout(devnull):
            genotypes = ProbabilisticGenotypes(genotype_names=samps)
            genotypes.add_vcf(vcf)

            snps = count_snps(bamfile_location=bam, chromosome2positions=genotypes.get_chromosome2positions(),
                barcode_handler=barcode_handler, joblib_n_jobs=1)

            likelihoods, posterior_probabilities = Demultiplexer.predict_posteriors(
                snps, genotypes=genotypes, barcode_handler=barcode_handler, doublet_prior=0.5)

    refmt_df(posterior_probabilities)
   
if __name__ == '__main__':
    sys.exit(main(sys.argv))
