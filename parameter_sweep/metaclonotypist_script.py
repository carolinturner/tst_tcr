import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


import igraph
import seaborn as sns
import scipy.stats

import pyrepseq as prs
import pyrepseq.plotting as pp
from metaclonotypist import *

plt.style.use('seaborn-v0_8-paper')

print(snakemake.wildcards)
chain = snakemake.wildcards['chain']
mincount = int(snakemake.wildcards['mincount'])
max_edits = 2
max_tcrdist = int(snakemake.wildcards['max_tcrdist'])
clustering = snakemake.wildcards['clustering']
if clustering == 'multilevel':
    clustering_kwargs = dict(resolution=10)
elif clustering == 'leiden':
    clustering_kwargs = dict(resolution=0.1,
                         objective_function='CPM',
                         n_iterations=4)
min_donors = 4
testmethod = snakemake.wildcards['testmethod']

chain_letter = chain[0].upper()
df = pd.read_csv(f'/Users/rishikasaxena/TRACERx_TCR_signatures-1/TST/combined_subsampled_5000_10000_{chain}_withduplicates_processed.csv')

hlas = pd.read_csv('/Users/rishikasaxena/TRACERx_TCR_signatures-1/TST/hlas.csv', index_col=0)
hlas = flatten_hlas(hlas)


# filter clones < mincount
df = df[df['clonal_count']>=mincount]
# only keep samples found in both datasets
df = df[df['Sample.ID'].isin(hlas.index)]
hlas = hlas.loc[list(set(df['Sample.ID']))]
# filter hlas < min_donors
hlas = hlas[hlas.columns[hlas.sum(axis=0)>=min_donors]]
# filter class I HLAS
hlas = hlas[hlas.columns[hlas.columns.str.startswith('D')]]


clusters = metaclonotypist(df, chain=chain,
                                   max_tcrdist=max_tcrdist, max_edits=max_edits,
                                   clustering=clustering, clustering_kwargs=clustering_kwargs)
clusters['Sample.ID'] = df.loc[clusters.index]['Sample.ID']


ndonors = clusters.groupby('cluster').apply(lambda cluster: len(cluster['Sample.ID'].unique()))
clusters = clusters[clusters['cluster'].isin(ndonors[(ndonors >= min_donors)].index)]
   
cluster_association = hla_association(clusters, hlas, method=testmethod)


hla_metaclones = cluster_association[cluster_association['significant']]

params_str = f'{chain}_td{max_tcrdist}_mc{mincount}_{testmethod}_{clustering}'
#clusters.to_csv(f'data/metaclonotypist/clusters_{params_str}.csv')
#cluster_association.to_csv(f'data/metaclonotypist/clusters_association_{params_str}.csv', index=False)
#hla_metaclones.to_csv(f'data/metaclonotypist/clusters_association_significant_{params_str}.csv', index=False)

sig_clusters = clusters[(clusters['cluster'].isin(hla_metaclones['cluster']))].reset_index()
sig_clusters = sig_clusters.merge(hla_metaclones, on='cluster')
hla_match = [hlas.loc[row['Sample.ID']][row['hla']] for ind, row in sig_clusters.iterrows()]
sig_clusters = sig_clusters.iloc[hla_match]
# sig_clusters.to_csv(f'data/metaclonotypist/sig_clusters_{params_str}.csv')


# shuffle hlas
hlas_shuffled = hlas.copy()
hlas_shuffled.index = np.random.permutation(hlas_shuffled.index)
cluster_association_shuffled = hla_association(clusters, hlas_shuffled, method=testmethod)
#cluster_association_shuffled.to_csv(f'data/metaclonotypist/clusters_association_shuffled_{params_str}.csv', index=False)
hla_metaclones_shuffled = cluster_association_shuffled[cluster_association_shuffled['significant']]
sig_clusters_shuffled = clusters[(clusters['cluster'].isin(hla_metaclones_shuffled['cluster']))].reset_index()
sig_clusters_shuffled = sig_clusters_shuffled.merge(hla_metaclones_shuffled[['cluster', 'hla']],on='cluster',how='left')
hla_match_shuffled = [hlas.loc[row['Sample.ID']][row['hla']] for ind, row in sig_clusters_shuffled.iterrows()]
sig_clusters_shuffled = sig_clusters_shuffled.iloc[hla_match_shuffled]

data = [['chain', chain],
        ['max_tcrdist', max_tcrdist],
        ['mincount', mincount],
        ['testmethod', testmethod],
        ['clustering', clustering],
        ['nassociations', len(hla_metaclones)],
        ['nassociations_shuffled', len(hla_metaclones_shuffled)],
        ['nmetaclones', len(hla_metaclones['cluster'].unique())],
        ['nmetaclones_shuffled', len(hla_metaclones_shuffled['cluster'].unique())],
        ['clustered_fraction', len(clusters)/len(df)],
        ['sig_clonotype_fraction', len(sig_clusters)/len(df)],
        ['sig_clonotype_fraction_shuffled', len(sig_clusters_shuffled)/len(df)],
        ['sig_read_fraction', df.loc[sig_clusters['index']]['clonal_count'].sum()/df['clonal_count'].sum()],
        ['sig_read_fraction_shuffled', df.loc[sig_clusters_shuffled['index']]['clonal_count'].sum()/df['clonal_count'].sum()],
        ['id_fraction', len(sig_clusters['Sample.ID'].unique())/len(df['Sample.ID'].unique())],
        ['id_fraction_shuffled', len(sig_clusters_shuffled['Sample.ID'].unique())/len(df['Sample.ID'].unique())],
       ]
index, values = list(zip(*data))
pd.Series(index=index, data=values, name='results').to_csv(f'/Users/rishikasaxena/TRACERx_TCR_signatures-1/TST/results/stats_{params_str}.csv', index=True)

