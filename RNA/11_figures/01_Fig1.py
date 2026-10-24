#!/usr/bin/env python3
"""Plot the RNA atlas, marker expression and sample composition."""

import json
import scanpy as sc
import pandas as pd
from plot_utils import proj_dir, palette, save, umap, dotplot

data = sc.read_h5ad(proj_dir / 'data/RNA/output/annotation/human_final.h5ad')
definition = json.loads((proj_dir / 'RNA/00_reference/cell_types.json').read_text())
umap(data,'cell_type','Fig1b_RNA_atlas')
genes = [g for ct in definition['major'] for g in definition['markers'][ct]]
dotplot(data,genes,'cell_type','Fig1d_RNA_markers')
for group in ['sample','disease']:
    table = pd.crosstab(data.obs[group],data.obs.cell_type,normalize='index')
    ax = table.plot.bar(stacked=True,color=[palette.get(x,'#999999') for x in table.columns],figsize=(6,3))
    ax.set(ylabel='Cell fraction',xlabel='')
    ax.legend(bbox_to_anchor=(1,1),frameon=False)
    save(ax.figure,f'SFig1_{group}_composition')
for gene in ['OTX2','RGS5','ABCA9','KCNK2','SLC4A4','FLT1','CYBB','CFAP299','GFAP','UGT8']:
    if gene in data.var_names:
        fig = sc.pl.umap(data,color=gene,show=False,return_fig=True)
        save(fig,f'SFig1_marker_{gene}')
