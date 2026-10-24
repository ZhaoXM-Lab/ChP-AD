#!/usr/bin/env python3
"""Plot cCRE burdens, NMF modules, epigenomic scores and functional concordance."""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
base = proj_dir / 'data/ATAC/output/ccres'
out = proj_dir / 'figures/ATAC'
out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':7,'pdf.fonttype':42})
def save(fig,name):
    fig.savefig(out / f'{name}.pdf',bbox_inches='tight')
    plt.close(fig)
diff = pd.read_csv(base / 'differential_ccres.csv')
for name,ct in [('major',~diff.cell_type.isin(['CPEC-1','CPEC-2'])),('CPEC',diff.cell_type.isin(['CPEC-1','CPEC-2']))]:
    selected = diff.loc[ct & diff.cell_type_specific].drop_duplicates(['cell_type','direction','peak']).copy()
    selected['specific'] = ~selected.duplicated(['peak','direction'],keep=False)
    table = selected.groupby(['cell_type','direction','specific']).size().unstack(['direction','specific'],fill_value=0)
    fig,ax = plt.subplots(figsize=(6,3))
    bottom = pd.Series(0,index=table.index)
    for direction,color in [('AD','#B02425'),('WT','#5355A2')]:
        for specific in [True,False]:
            values = table.get((direction,specific),pd.Series(0,index=table.index))
            ax.bar(table.index,values,bottom=bottom,color=color if specific else '#bbbbbb',label=f'{direction}: '+('specific' if specific else 'shared'))
            bottom += values
    ax.tick_params(axis='x',rotation=90)
    ax.set_ylabel('Differential cell-specific cCREs')
    ax.legend(frameon=False)
    save(fig,f'Fig3c_{name}')
h = pd.read_csv(base / 'nmf/normalized_activity.tsv',sep='\t',index_col=0)
fig,ax = plt.subplots(figsize=(8,3))
sns.heatmap(h,cmap='Reds',ax=ax)
save(fig,'Fig3d_NMF_activity')
score = pd.read_csv(base / 'epigenomic_information.csv')
fig,ax = plt.subplots(figsize=(8,3))
sns.boxplot(data=score,x='cell_type',y='score',hue='disease',showfliers=False,ax=ax)
ax.tick_params(axis='x',rotation=90)
save(fig,'Fig3d_epigenomic_information')
module = pd.read_csv(base / 'nmf/module_cell_accessibility.csv')
for m,z in module.groupby('module'):
    fig,ax = plt.subplots(figsize=(7,3))
    sns.boxplot(data=z,x='cell_type',y='accessibility',hue='disease',showfliers=False,ax=ax)
    ax.tick_params(axis='x',rotation=90)
    save(fig,f'Fig3d_{m}_accessibility')
terms = pd.read_csv(base / 'functional_concordance.csv')
fig,ax = plt.subplots(figsize=(7,3))
sns.barplot(data=terms,x='cell_type',y='overlap_ratio',hue='disease',ax=ax)
ax.tick_params(axis='x',rotation=90)
save(fig,'Fig3e_functional_concordance')
trends = pd.read_csv(base / 'nmf/rank_diagnostics.csv')
fig,axes = plt.subplots(1,2,figsize=(6,3))
for ax,col in zip(axes,['mean_sparsity','mean_entropy']):
    ax.plot(trends.k,trends[col],marker='o'); ax.set(xlabel='NMF rank',ylabel=col)
save(fig,'NMF_rank_diagnostics')
for file in sorted((base / 'great').glob('*.csv')):
    table = pd.read_csv(file).sort_values('p_adjust_hyper').head(10)
    if len(table):
        ax = table.sort_values('fold_enrichment').plot.barh(x='description',y='fold_enrichment',legend=False,figsize=(6,3))
        save(ax.figure,'GREAT_'+file.stem)
