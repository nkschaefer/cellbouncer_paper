#! /usr/bin/env python3
import sys
import os
from collections import Counter, defaultdict
import numpy as np
import pandas as pd
from pomegranate import *
import scipy
from betabinomial import BetaBinomialDistribution
import math

def fit_model_2part(sub, init_phi=10):
    """
    This function is for when there are only 2 models -- both mito or one mito.
    This is appropriate for within-species fusions.
    """

    # sub1 = count1 > count2 (het or ID1)
    # sub3 = flip around count 1 < count 2
    sub1 = sub.loc[sub['count1'] >= sub['count2'],:]
    sub2 = sub.loc[sub['count2'] > sub['count1'],:]
    sub3 = pd.DataFrame({'bc': sub2['bc'], 'id1': sub2['id2'], 'id2': sub2['id1'],
        'count1': sub2['count2'], 'count2': sub2['count1'],
        'lib': sub2['lib'], 'species': sub2['species'], 'tot': sub2['tot']})
    sub = pd.concat([sub1, sub3], axis=0)
    sub['name_het'] = np.where((sub['id1'].astype('str') < sub['id2'].astype('str')), sub['id1'] + "+" + sub['id2'], 
        sub['id2'] + "+" + sub['id1'])
   
    m = np.mean(sub['count1']/sub['tot'])

    dists = [BetaBinomialDistribution(m*init_phi, (1.0-m)*init_phi), \
            BetaBinomialDistribution(0.99*init_phi, 0.01*init_phi)]

    dat = np.hstack([sub['count1'].values.reshape(-1,1), sub['tot'].values.reshape(-1,1)])

    mod = GeneralMixtureModel(dists)
    mod = mod.fit(dat)
    
    probs = mod.predict_log_proba(dat)
    cols = {'bc': sub['bc']}

    a0 = mod.distributions[0].alpha
    b0 = mod.distributions[0].beta
    a1 = mod.distributions[1].alpha
    b1 = mod.distributions[1].beta
   
    df = pd.DataFrame({'bc': sub['bc'], 'het': probs[:,0], 'hom': probs[:,1]})
    if a0/(a0+b0) > 0.95 and a1/(a1+b1) > 0.95:
        print("BOTH HOM", file=sys.stderr)
        df['het'] = math.log(1e-4)
        df['hom'] = math.log(1.0-1e-4)
    elif abs(a0/(a0+b0) - a1/(a1+b1)) < 0.1:
        print("BOTH HET", file=sys.stderr)
        df['het'] = math.log(1.0-1e-4)
        df['hom'] = math.log(1e-4)
    elif a0/(a0+b0) > a1/(a1+b1):
        print("SWITCH", file=sys.stderr)
        df['het'] = probs[:,1]
        df['hom'] = probs[:,0]
    df['name_het'] = sub['name_het']
    df['name_hom'] = sub['id1']
    
    df['name'] = ''
    df.loc[df['het'] > df['hom'],'name'] = df.loc[df['het'] > df['hom'],'name_het']
    df.loc[df['hom'] > df['het'],'name'] = df.loc[df['hom'] > df['het'],'name_hom']
    
    return df.drop(['name_het', 'name_hom'], axis=1)

def fit_model(sub, names, pre_fit_mod=None, init_phi=10, excl_low=False):
    """
    This applies for when there are 3 models - species1, species2, or both
    This is appropriate for inter-species fusions.
    """

    name1, name2 = names
    hetname = '{}+{}'.format(sorted(names)[0], sorted(names)[1])
    
    dists = [BetaBinomialDistribution(0.01*init_phi, 0.99*init_phi), \
            BetaBinomialDistribution(0.5*init_phi, 0.5*init_phi), \
            BetaBinomialDistribution(0.99*init_phi, 0.01*init_phi)]
    dists[0].name = name1
    dists[1].name = hetname
    dists[2].name = name2
    
    if excl_low:
        del dists[0]

    dat = np.hstack([sub['count1'].values.reshape(-1,1), sub['tot'].values.reshape(-1,1)])
    
    mod = None
    if pre_fit_mod is not None:
        mod = pre_fit_mod
    else:
        mod = GeneralMixtureModel(dists)
        mod = mod.fit(dat)
    
    print(mod)
    probs = mod.predict_log_proba(dat)
    probs2 = mod.predict_proba(dat)

    cols = {'bc': sub['bc']}
    cols2 = {'bc': sub['bc'], 'gene': sub['gene']}
    
    a0 = mod.distributions[0].alpha
    b0 = mod.distributions[0].beta
    a1 = mod.distributions[1].alpha
    b1 = mod.distributions[1].beta
    a2 = mod.distributions[2].alpha
    b2 = mod.distributions[2].beta

    if excl_low:
        cols[hetname] = probs[:,0]
        cols2[hetname] = probs2[:,0]
        cols[name2] = probs[:,1]
        cols2[name2] = probs2[:,1]
    else:
        cols[name1] = probs[:,0]
        cols2[name1] = probs2[:,0]
    
        cols[hetname] = probs[:,1]
        cols2[hetname] = probs2[:,1]

        cols[name2] = probs[:,2]
        cols2[name2] = probs2[:,2]
        
    probdf = pd.DataFrame(cols)
    probdf2 = pd.DataFrame(cols2)
    
    probdf2['id'] = probdf2.drop(['bc', 'gene'], axis=1).idxmax(axis=1)
    probdf2['p'] = probdf2.drop(['bc', 'gene', 'id'], axis=1).max(axis=1)

    prob_cell = probdf.groupby('bc').sum()
    prob_cell['id'] = prob_cell.idxmax(axis=1)
    
    prob_cell['max'] = prob_cell.drop(['id'], axis=1).max(axis=1)
    
    if excl_low:
        prob_cell['p'] = 1.0 / (np.exp(prob_cell[hetname] - prob_cell['max']) + \
            np.exp(prob_cell[name2] - prob_cell['max']))
    else:
        prob_cell['p'] = 1.0 / (np.exp(prob_cell[name1] - prob_cell['max']) + \
            np.exp(prob_cell[hetname] - prob_cell['max']) + \
            np.exp(prob_cell[name2] - prob_cell['max']))
    
    return (mod, prob_cell, probdf2)
    

def main(args):
    
    # Intra species
    tabdf = pd.read_csv('intra_species_dat.tsv', sep='\t')
    
    tabdf['tot'] = tabdf['count1'] + tabdf['count2']
    
    cell_probs_H_nofix = []
    gene_probs_H_nofix = []
    cell_probs_C_nofix = []
    gene_probs_C_nofix = []
    
    subH = tabdf.loc[(tabdf['species'] == 'Human_Human'),:]
    subC = tabdf.loc[(tabdf['species'] == 'Chimp_Chimp'),:]
    
    subHC1 = subH.groupby(['bc', 'id1', 'id2', 'lib', 'species'])['count1'].sum().reset_index()
    subHC2 = subH.groupby(['bc'])['count2'].sum().reset_index()
    subH = subHC1.merge(subHC2, left_on='bc', right_on='bc')

    subCC1 = subC.groupby(['bc', 'id1', 'id2', 'lib', 'species'])['count1'].sum().reset_index()
    subCC2 = subC.groupby(['bc'])['count2'].sum().reset_index()
    subC = subCC1.merge(subCC2, left_on='bc', right_on='bc')
    subH['tot'] = subH['count1'] + subH['count2']
    subC['tot'] = subC['count1'] + subC['count2']
    
    res_H_cell = fit_model_2part(subH)
    res_C_cell = fit_model_2part(subC)
    
    cell_probs_H_nofix.append(res_H_cell)
    cell_probs_C_nofix.append(res_C_cell)
        
    cell_probs_H_nofix_df = pd.concat(cell_probs_H_nofix, axis=0)
    cell_probs_C_nofix_df = pd.concat(cell_probs_C_nofix, axis=0)

    cell_probs_H_nofix_df['species'] = 'Human_Human'
    cell_probs_C_nofix_df['species'] = 'Chimp_Chimp'

    cell_probs_nofix_df = pd.concat([cell_probs_H_nofix_df, cell_probs_C_nofix_df], axis=0)
    
    cell_probs_nofix_df.index = cell_probs_nofix_df['bc']
    cell_probs_nofix_df.drop(['bc'], axis=1, inplace=True)
    cell_probs_nofix_df.index.name = 'bc'
    cell_probs_nofix_df.to_csv('mito_ase_results_intra_species.tsv', sep='\t')

    cutoff = 0.5
    
    # Inter species
    tab = pd.read_csv('inter_species_dat.tsv', sep='\t')
    
    tab['frac'] = (tab['count1'])/(tab['count1'] + 1 + tab['count2'])
    tab['tot'] = tab['count1'] + tab['count2']
    
    # Chimp/Bonobo
    subcb = tab.loc[(tab['s1'] == 'Chimp') & \
        (tab['s2'] == 'Bonobo') & \
        (tab['species'] == 'Chimp_Bonobo'),:]
    
    subcb_C = tab.loc[(tab['s1'] == 'Chimp') & \
        (tab['s2'] == 'Bonobo') & \
        (tab['species'] == 'Chimp_Chimp'),:]
    
    modcb, cell_probs_cb, gene_probs_cb = fit_model(subcb, ['Bonobo', 'Chimp'], excl_low=False)
    modcb2, cell_probs_cb_c, gene_probs_cb_c = fit_model(subcb_C, ['Bonobo', 'Chimp'], \
        pre_fit_mod=modcb)
    
    print("Inter-species C/B")
    print(cell_probs_cb.loc[cell_probs_cb['p'] > cutoff,'id'].value_counts() /\
        cell_probs_cb.loc[cell_probs_cb['p'] > cutoff,:].shape[0])
    print("Inter-species C/B, Chimp cells")
    print(cell_probs_cb_c.loc[cell_probs_cb_c['p'] > cutoff,'id'].value_counts() /\
        cell_probs_cb_c.loc[cell_probs_cb_c['p'] > cutoff,:].shape[0])
   
    
    cell_probs_cb.loc[cell_probs_cb['p'] > cutoff,:].to_csv("mito_ase_results_CB.tsv", sep="\t")
    gene_probs_cb.loc[gene_probs_cb['p'] > cutoff,:].to_csv('mito_ase_results_CB_gene.tsv', sep='\t')

    # Human/Chimp

    init_phi = 10
    
    cell_probs_HC_all = []
    gene_probs_HC_all = []

    cell_probs_H_all = []
    gene_probs_H_all = []

    cell_probs_C_all = []
    gene_probs_C_all = []


    sub = tab.loc[(tab['s1'] == 'Human') & \
            (tab['s2'] == 'Chimp') & \
            (tab['species'] == 'Human_Chimp'),:]
    
    sub_C = tab.loc[(tab['s1'] == 'Human') & \
            (tab['s2'] == 'Chimp') & \
            (tab['species'] == 'Chimp_Chimp'),:]

    sub_H = tab.loc[(tab['s1'] == 'Human') & \
            (tab['s2'] == 'Chimp') & \
            (tab['species'] == 'Human_Human'),:]
    mod, cell_probs, gene_probs = fit_model(sub, ['Chimp', 'Human'])
    
    mod2, cell_probs_C, gene_probs_C = fit_model(sub_C, ['Chimp', 'Human'], \
        pre_fit_mod=mod)
    
    mod3, cell_probs_H, gene_probs_H = fit_model(sub_H, ['Chimp', 'Human'], \
        pre_fit_mod=mod)
    
    cell_probs_HC_all.append(cell_probs)
    gene_probs_HC_all.append(gene_probs)
    cell_probs_H_all.append(cell_probs_H)
    gene_probs_H_all.append(gene_probs_H)
    cell_probs_C_all.append(cell_probs_C)
    gene_probs_C_all.append(gene_probs_C)

    cell_probs_HC_all_df = pd.concat(cell_probs_HC_all, axis=0)
    gene_probs_HC_all_df = pd.concat(gene_probs_HC_all, axis=0)
    cell_probs_H_all_df = pd.concat(cell_probs_H_all, axis=0)
    gene_probs_H_all_df = pd.concat(gene_probs_H_all, axis=0)
    cell_probs_C_all_df = pd.concat(cell_probs_C_all, axis=0)
    gene_probs_C_all_df = pd.concat(gene_probs_C_all, axis=0)
    
    props = cell_probs_HC_all_df.loc[cell_probs_HC_all_df['p'] > cutoff,'id'].value_counts() /\
        cell_probs_HC_all_df.loc[cell_probs_HC_all_df['p'] > cutoff,:].shape[0]

    props_H = cell_probs_H_all_df.loc[cell_probs_H_all_df['p'] > cutoff,'id'].value_counts() /\
        cell_probs_H_all_df.loc[cell_probs_H_all_df['p'] > cutoff,:].shape[0]

    props_C = cell_probs_C_all_df.loc[cell_probs_C_all_df['p'] > cutoff, 'id'].value_counts() /\
        cell_probs_C_all_df.loc[cell_probs_C_all_df['p'] > cutoff,:].shape[0]
    
    cell_probs_HC_all_df.loc[cell_probs_HC_all_df['p'] > cutoff,:].to_csv('mito_ase_results_HC.tsv', sep='\t')
    gene_probs_HC_all_df.loc[gene_probs_HC_all_df['p'] > cutoff,:].to_csv('mito_ase_results_HC_gene.tsv', sep='\t')
    
    print("Inter-species H/C")
    print(props)
    print("Inter-species H/C Human")
    print(props_H)
    print("Inter-species H/C Chimp")
    print(props_C)
    
    print(cell_probs_HC_all_df)


if __name__ == '__main__':
    sys.exit(main(sys.argv))
