#! /usr/bin/env python3
import sys
import os
import anndata as ad
import pandas as pd
from solo import hashsolo

def main(args):
    if len(args) < 2:
        print("USAGE: run_hashsolo.py <MTX_DIR>", file=sys.stderr)
        exit(1)

    adata = ad.read_mtx('{}/matrix.mtx.gz'.format(args[1]))
    adata = adata.transpose()
    features = pd.read_csv('{}/features.tsv.gz'.format(args[1]), names=[0])
    barcodes = pd.read_csv('{}/barcodes.tsv.gz'.format(args[1]), names=[0])
    
    adata.obs.index = barcodes[0].values
    adata.var.index = features[0].values
    
    #adata.write_h5ad('TEST.h5ad')

    hashsolo.hashsolo(adata)
    
    df = adata.obs.drop(['most_likely_hypothesis', 'cluster_feature', 'negative_hypothesis_probability', 'singlet_hypothesis_probability', 'doublet_hypothesis_probability'], axis=1)
    df.loc[df['Classification'] == 'Doublet','Classification'] = 'multiplet'
    df = df.loc[df['Classification'] != 'Negative',:]    
    
    df.to_csv(sys.stdout, sep='\t', header=False)



if __name__ == '__main__':
    sys.exit(main(sys.argv))
