#!/usr/bin/env python3
"""Shared figure styling and marker dot plots."""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'figures/RNA'
out.mkdir(parents=True, exist_ok=True)
plt.rcParams.update({'font.size':7,'pdf.fonttype':42,'ps.fonttype':42,'savefig.dpi':300})
palette = {'CPEC':'#BC80BD','CPEC-1':'#B174E7','CPEC-2':'#BEA1F7','Ependymal':'#637939',
           'Astrocyte':'#A6D854','Oligodendrocyte':'#CEDB9C','Macrophage':'#5254A3','Endothelial':'#FFED6F',
           'Pericyte':'#AD494A','V-SMC':'#E7969C','M-SMC':'#843C39','Fibroblast':'#D6616B',
           'OPC':'#D9F4CD','Neuron':'#FFD21F','CD8+T':'#FF2828'}

def save(fig, name):
    fig.savefig(out / f'{name}.pdf', bbox_inches='tight')
    plt.close(fig)

def umap(data, column, name):
    fig, ax = plt.subplots(figsize=(4,4))
    for label in sorted(data.obs[column].unique()):
        mask = (data.obs[column] == label).to_numpy()
        xy = data.obsm['X_umap'][mask]
        ax.scatter(xy[:,0],xy[:,1],s=0.5,label=label,c=palette.get(label),rasterized=True)
    ax.set(xlabel='UMAP1',ylabel='UMAP2',xticks=[],yticks=[])
    ax.legend(markerscale=5,fontsize=6,bbox_to_anchor=(1,1),frameon=False)
    save(fig,name)

def dotplot(data, genes, column, name):
    genes = [g for g in genes if g in data.var_names]
    groups = list(data.obs[column].cat.categories) if isinstance(data.obs[column].dtype, pd.CategoricalDtype) else sorted(data.obs[column].unique())
    avg, pct = [], []
    for group in groups:
        x = data[data.obs[column] == group,genes].X
        avg.append(np.asarray(x.mean(axis=0)).ravel())
        pct.append(np.asarray((x > 0).mean(axis=0)).ravel())
    avg = np.array(avg)
    z = (avg-avg.mean(axis=0))/np.maximum(avg.std(axis=0),1e-12)
    xx,yy = np.meshgrid(np.arange(len(genes)),np.arange(len(groups)))
    fig,ax = plt.subplots(figsize=(max(4,len(genes)*0.28),max(2,len(groups)*0.25)))
    points = ax.scatter(xx.ravel(),yy.ravel(),s=np.array(pct).ravel()*45,c=z.ravel(),cmap='RdBu_r',vmin=-2,vmax=2)
    ax.set(xticks=range(len(genes)),xticklabels=genes,yticks=range(len(groups)),yticklabels=groups)
    ax.tick_params(axis='x',rotation=90)
    fig.colorbar(points,ax=ax,label='Mean expression z-score')
    save(fig,name)
