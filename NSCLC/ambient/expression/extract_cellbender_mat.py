#! /usr/bin/env python3
import sys
import os
import scipy
import gzip
from cellbender.remove_background.downstream import *

def main(args):
    if len(args) < 3:
        print("Please provide file.h5 outdir", file=sys.stderr)
        exit(1)
    
    inh5 = args[1]
    outdir = args[2]
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    a = anndata_from_h5(inh5)
    mat = a.X.transpose()
    bcs = a.obs.index
    f = open('{}/barcodes.tsv'.format(outdir), 'w')
    for bc in bcs:
        print(bc, file=f)
    f.close()

    f = open('{}/features.tsv'.format(outdir), 'w')
    for feat in a.var.index:
        print(feat, file=f)
    f.close()

    scipy.io.mmwrite("{}/matrix.mtx".format(outdir), mat)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
