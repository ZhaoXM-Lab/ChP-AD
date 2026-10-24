#!/usr/bin/env python3
"""Plot hormone differences and prioritized drug-gene interactions."""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plot_utils import proj_dir, save

base = proj_dir / 'data/RNA/output/hormone'
for role in ['HPC','HRC']:
    diff = pd.read_csv(base / f'{role}_difference.csv')
    matrix = diff.pivot_table(index='Hormone',columns='cell_type',values='delta')
    fig,ax = plt.subplots(figsize=(7,max(3,len(matrix)*0.15)))
    sns.heatmap(matrix,cmap='RdBu_r',center=0,ax=ax)
    save(fig,f'SFig7_{role}_difference')
edges = pd.read_csv(base / 'network_difference.csv')
edges['pair'] = edges.sender + ' -> ' + edges.receiver
matrix = edges.pivot_table(index='Hormone',columns='pair',values='delta')
fig,ax = plt.subplots(figsize=(max(8,len(matrix.columns)*0.14),max(3,len(matrix)*0.2)))
sns.heatmap(matrix,cmap='RdBu_r',center=0,ax=ax)
save(fig,'Fig4g_hormone_network')
for mode, group in edges.groupby('mode'):
    matrix = group.pivot_table(index='Hormone',columns='pair',values='delta')
    fig,ax = plt.subplots(figsize=(max(6,len(matrix.columns)*0.14),max(3,len(matrix)*0.2)))
    sns.heatmap(matrix,cmap='RdBu_r',center=0,ax=ax)
    save(fig,f'SFig7c_{mode}')
drug = pd.read_csv(proj_dir / 'data/RNA/output/drug2cell/drug_gene_interactions.csv')
for name, selected in [('SFig11',drug),('Fig6f',drug.loc[drug.cell_type.isin(['CPEC','CPEC-1','CPEC-2','Macrophage'])])]:
    matrix = selected.pivot_table(index='drug',columns=['cell_type','gene'],values='score')
    if matrix.empty:
        continue
    fig,ax = plt.subplots(figsize=(max(6,len(matrix.columns)*0.18),max(3,len(matrix)*0.13)))
    sns.heatmap(matrix,cmap='Reds',ax=ax)
    save(fig,f'{name}_drug_gene')
