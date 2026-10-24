#!/usr/bin/env python3
"""Plot human and mouse CPEC states, functional markers and enrichment."""

import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
from plot_utils import proj_dir, save, umap, dotplot, palette

for species in ['human','mouse']:
    data = sc.read_h5ad(proj_dir / f'data/RNA/output/annotation/{species}_cpec.h5ad')
    umap(data,'cell_type',f'{species}_CPEC_umap')
    if species == 'human':
        genes = ['DNAH5','SLC24A4','NDUFA1','NDUFAB1']
    else:
        homologues = pd.read_csv(proj_dir / 'data/RNA/output/annotation/human_mouse_orthologs.tsv',sep='\t')
        all_genes = homologues.mmusculus_homolog_associated_gene_name.dropna().unique().tolist()
        genes = homologues.loc[homologues.external_gene_name.isin(['DNAH5','SLC24A4','NDUFA1','NDUFAB1']), 'mmusculus_homolog_associated_gene_name'].dropna().unique().tolist()
        atlas = sc.read_h5ad(proj_dir / 'data/RNA/output/annotation/mouse_final.h5ad')
        umap(atlas,'cell_type','SFig3c_mouse_atlas')
        dotplot(atlas,all_genes,'cell_type','SFig3d_mouse_markers')
    dotplot(data,genes,'cell_type',f'{species}_CPEC_marker_dotplot')
    for gene in genes:
        if gene in data.var_names:
            fig = sc.pl.umap(data,color=gene,show=False,return_fig=True)
            save(fig,f'{species}_CPEC_feature_{gene}')
    composition = pd.crosstab(data.obs['sample'],data.obs.cell_type,normalize='index')
    ax = composition.plot.bar(stacked=True,color=[palette[x] for x in composition.columns],figsize=(5,3))
    save(ax.figure,f'{species}_CPEC_sample_composition')
markers = pd.read_csv(proj_dir / 'data/RNA/output/annotation/human_cpec_markers.csv')
markers = markers.loc[markers.group == 'CPEC-1']
fig,ax = plt.subplots(figsize=(4,3))
sig = (markers.pvals_adj < 0.05) & (markers.logfoldchanges.abs() > 0.3)
import numpy as np
ax.scatter(markers.logfoldchanges,-np.log10(markers.pvals_adj.clip(lower=1e-300)),c=np.where(sig,'#B174E7','#bbbbbb'),s=2,rasterized=True)
ax.set(xlabel='log2 fold change',ylabel='-log10 adjusted P')
save(fig,'SFig3b_CPEC_volcano')
for ct in ['CPEC-1','CPEC-2']:
    base = proj_dir / 'data/RNA/output/enrichment'
    report = next((base / f'{ct}_gsea').glob('*prerank.report*.csv'))
    gsea = pd.read_csv(report).sort_values('NES',ascending=False).head(10).sort_values('NES')
    ax = gsea.plot.barh(x='Term',y='NES',legend=False,figsize=(6,3),color=palette[ct])
    save(ax.figure,f'Fig2_{ct}_GSEA')
    enrich = pd.read_csv(base / f'{ct}_markers/GO_Biological_Process_2025.human.enrichr.reports.txt',sep='\t')
    enrich = enrich.loc[enrich['Adjusted P-value'] < 0.05].sort_values('Combined Score').tail(10)
    ax = enrich.plot.barh(x='Term',y='Combined Score',legend=False,figsize=(6,3),color=palette[ct])
    save(ax.figure,f'SFig4a_{ct}_Enrichr')
