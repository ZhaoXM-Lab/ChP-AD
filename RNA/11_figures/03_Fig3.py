#!/usr/bin/env python3
"""Plot consensus DEG burdens and expression changes of specific TFs."""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from plot_utils import proj_dir, save

deg = pd.read_csv(proj_dir / 'data/RNA/output/differential_expression/consensus_degs.csv')
for name, labels in [('major',~deg.cell_type.isin(['CPEC-1','CPEC-2'])),('CPEC',deg.cell_type.isin(['CPEC-1','CPEC-2']))]:
    counts = deg.loc[labels].groupby(['cell_type','direction','specific']).size().unstack(['direction','specific'],fill_value=0)
    fig,ax = plt.subplots(figsize=(7,3))
    bottom = pd.Series(0,index=counts.index)
    for direction,color in [('Up','#B02425'),('Down','#5355A2')]:
        for specific in [True,False]:
            values = counts.get((direction,specific),pd.Series(0,index=counts.index))
            ax.bar(counts.index,values,bottom=bottom,color=color if specific else '#bbbbbb',label=f'{direction}: '+('specific' if specific else 'shared'))
            bottom += values
    ax.tick_params(axis='x',rotation=90)
    ax.set_ylabel('Consensus DEGs')
    ax.legend(frameon=False)
    save(fig,f'Fig3a_{name}_DEGs')
tf = pd.read_csv(proj_dir / 'data/ATAC/output/grn/specific_tfs.csv')
selected = deg.merge(tf[['cell_type','TF']].drop_duplicates(),left_on=['cell_type','gene'],right_on=['cell_type','TF'])
if len(selected):
    matrix = selected.pivot_table(index='gene',columns='cell_type',values='avg_log2FC')
    fig,ax = plt.subplots(figsize=(6,max(3,len(matrix)*0.12)))
    sns.heatmap(matrix,cmap='RdBu_r',center=0,ax=ax,cbar_kws={'label':'AD-control log2 fold change'})
    save(fig,'Fig3b_differential_TFs')
for file in sorted((proj_dir / 'data/RNA/output/enrichment').glob('*_*/GO_Biological_Process_2025.human.enrichr.reports.txt')):
    if not file.parent.name.endswith(('_Up','_Down')):
        continue
    table = pd.read_csv(file,sep='\t')
    table = table.loc[table['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value').head(10)
    if table.empty:
        continue
    fig,ax = plt.subplots(figsize=(6,3))
    import numpy as np
    points = ax.scatter(-np.log10(table['Adjusted P-value'].clip(lower=1e-300)),table.Term,
                        s=table['Odds Ratio']*3,c=table['Combined Score'],cmap='Reds')
    ax.set_xlabel('-log10 adjusted P')
    fig.colorbar(points,ax=ax,label='Combined score')
    save(fig,'SFig5_RNA_'+file.parent.name)
