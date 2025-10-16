#! /usr/bin/env python3
import sys
import os
from pomegranate import *
import pandas as pd
import numpy as np
import time

def PoissonGaussian(sgRNA: str, df, UMI_threshold):
    start = time.time()
    print(sgRNA)
    cell_barcodes = df[sgRNA].where(lambda x: x > UMI_threshold).dropna().index
    UMIs = df[sgRNA].where(lambda x: x > UMI_threshold).dropna().values
    log2_UMIs = np.log2(df[sgRNA].where(lambda x: x > UMI_threshold).dropna()).values.reshape(-1,1)
    
    maxtime = 3

    i = 0
    #array of 1000 linearly spaced points between -2 and max(log2_UMIs)+2
    gmm_x = np.linspace(-2, max(log2_UMIs)+2, 1000)
    #run until variable i != 0
    while i == 0 and time.time()-start < maxtime:
        #create a model using the GeneralMixtureModel.from_samples function, which is a mixture model of two distributions: PoissonDistribution and NormalDistribution
        model = GeneralMixtureModel.from_samples([PoissonDistribution, NormalDistribution], 2, log2_UMIs)
        #check that the model converges (no NaN) and that the Poisson distribution is the lower component by using an if-else loop
        #if there is any NaN in the probability of the model, set i to 0 to keep the loop running
        if numpy.isnan(model.probability(gmm_x)).any():
            i = 0
        else:
            #check that Poisson distribution is the lower component
            if model.distributions[0].name == 'PoissonDistribution':
                #if the first component of the model is PoissonDistribution and has a smaller parameter than the second component, set i to 1 to exit the loop
                if model.distributions[0].parameters[0] < model.distributions[1].parameters[0]:
                    i = 1    
                else:
                    i = 0
            #if the first component has a bigger parameter than the second, set i to 0 to keep the loop running
            elif model.distributions[0].parameters[0] > model.distributions[1].parameters[0]:
                i = 0

    #probability of the model is greater than 0.5 (second column is probability of belonging to second component)
    #return model.predict_proba(log2_UMIs)[:,1] > 0.5
    end = time.time()
    if end-start >= maxtime:
        return (False, pd.DataFrame({'cell_barcode': [], 'sg_ID': [], "UMIs": [], "PoissonGaussian": [] }))
    else:
        print('Time:', end - start)
        return (True, pd.DataFrame({'cell_barcode': cell_barcodes, 'sg_ID': sgRNA, 'UMIs': UMIs, 'PoissonGaussian': model.predict_proba(log2_UMIs)[:,1]}))

def main(args):
    if len(args) < 2:
        print("USAGE: ./pg.py [outprefix].counts", file=sys.stderr)
        exit(1)
    
    counts = pd.read_csv(args[1], sep='\t', index_col=0)
    results = []
    for sgRNA in counts.columns:
        try:
            success, df = PoissonGaussian(sgRNA, counts, 10)
            if success:
                results.append(df)
            else:
                print("Skip {}; time".format(sgRNA), file=sys.stderr)

        except:
            print("Skip {}; error".format(sgRNA), file=sys.stderr)
    
    results_all = np.vstack(results)
    results_all = pd.DataFrame(results_all, columns=['barcode', 'sgRNA', 'UMIs', 'prob'])
    results_all.index = results_all['barcode']
    results_all.drop(['barcode'], axis=1, inplace=True)
    
    results1 = results_all.loc[results_all['prob'] > 0.5,:]
    results1b = results1.drop(['UMIs', 'prob'], axis=1).\
        groupby(['barcode'], as_index=True).agg(lambda x: ','.join(list(x)))
    results1count = results1.drop(['UMIs', 'prob'], axis=1).\
        groupby(['barcode'], as_index=True).agg(lambda x: len(list(x)))
    
    results1 = results1b.merge(results1count, left_index=True, right_index=True)
    results1.columns = ['sgRNA', 'num']

    # Add in unassigned
    bc_missing = set(counts.index).difference(results1.index)
    bc_missing = list(bc_missing)
    results1c = pd.DataFrame({'barcode': bc_missing, 'sgRNA': ['WT']*len(bc_missing), 'num': [0]*len(bc_missing)})
    results1c.index = results1c['barcode']
    results1c = results1c.drop(['barcode'], axis=1)
    
    results1 = pd.concat([results1, results1c])
    
    results1.to_csv('{}.pg.assignments'.format(args[1].split('.counts')[0]), sep='\t', header=False)  

    results_all['has'] = 0
    results_all.loc[results_all['prob'] > 0.5,'has'] = 1
    
    # Add in unassigned
    bc_missing2 = set(counts.index).difference(results_all.index)
    bc_missing2 = list(bc_missing2)
    if len(bc_missing2) > 0:
        sgrnafill = results1.loc[results1['sgRNA'] != 'WT','sgRNA'][1].split(',')[0]
        results_all2 = pd.DataFrame({'barcode': bc_missing2, 'sgRNA': [sgrnafill]*len(bc_missing2), \
                'UMIs': [0.0]*len(bc_missing2), 'prob': 0.0, 'has': [0]*len(bc_missing2)})
        results_all2.index = results_all2['barcode']
        results_all2.drop(['barcode'], inplace=True, axis=1)
        results_all = pd.concat([results_all, results_all2])

    results_reshape = results_all.pivot(columns=['sgRNA'], values=['has'])
    for coln in results_reshape.columns:
        results_reshape.loc[results_reshape[coln].isna(),coln] = 0
        results_reshape[coln] = results_reshape[coln].astype('int')
    
    outn = args[1].split('.counts')[0] + '.pg.table'
    results_reshape.columns = results_reshape.columns.droplevel(0)
    #results_reshape.index += "-" + args[1].split('.counts')[0]
    results_reshape.to_csv(outn, sep='\t')
    

if __name__ == '__main__':
    sys.exit(main(sys.argv))
