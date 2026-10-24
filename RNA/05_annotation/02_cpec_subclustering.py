#!/usr/bin/env python3
"""Recluster epithelial nuclei and assign CPEC subtypes."""

import os
from pathlib import Path
import pandas as pd
import scanpy as sc

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
species = os.environ.get('SPECIES', 'human')
resolution = float(os.environ.get('CPEC_RESOLUTION', '0.1'))
out_dir = proj_dir / 'data/RNA/output/annotation'
atlas = sc.read_h5ad(out_dir / f'{species}_annotated.h5ad')
data = atlas[atlas.obs.cell_type == 'CPEC'].copy()
data.X = data.layers['counts'].copy()
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
data.raw = data.copy()
sc.pp.highly_variable_genes(data, flavor='seurat', n_top_genes=2000)
hvg = data[:, data.var.highly_variable].copy()
sc.pp.scale(hvg, max_value=10)
sc.pp.pca(hvg, n_comps=30, random_state=0)
sc.external.pp.harmony_integrate(hvg, key='sample', max_iter_harmony=30)
data.obsm['X_pca_harmony'] = hvg.obsm['X_pca_harmony']
sc.pp.neighbors(data, use_rep='X_pca_harmony', n_neighbors=10, random_state=0)
sc.tl.umap(data, random_state=0)
sc.tl.leiden(data, resolution=resolution, key_added='cpec_cluster', random_state=0)
sc.tl.rank_genes_groups(data, groupby='cpec_cluster', method='wilcoxon', pts=True)
sc.get.rank_genes_groups_df(data, group=None).to_csv(out_dir / f'{species}_cpec_cluster_markers.csv', index=False)
pd.crosstab(data.obs.cpec_cluster, data.obs['sample']).to_csv(out_dir / f'{species}_cpec_sample_counts.csv')
data.write_h5ad(out_dir / f'{species}_cpec_clusters.h5ad')
anno = pd.read_csv(proj_dir / 'RNA/00_reference/cpec_annotations.tsv', sep='\t', dtype=str)
anno = anno.loc[anno.species == species].set_index('cluster', verify_integrity=True)
if not set(data.obs.cpec_cluster.astype(str)).issubset(anno.index):
    raise ValueError('CPEC annotations do not cover all observed subclusters.')
labels = data.obs.cpec_cluster.astype(str)
include = labels.map(anno.include) == '1'
excluded = data.obs_names[~include]
data = data[include].copy()
data.obs['cell_type'] = labels.loc[include].map(anno.cell_type).astype('category')
if not set(data.obs.cell_type).issubset({'CPEC-1','CPEC-2'}):
    raise ValueError('Retained epithelial labels must be CPEC-1 or CPEC-2.')
sc.tl.rank_genes_groups(data, groupby='cell_type', method='wilcoxon', pts=True)
sc.get.rank_genes_groups_df(data, group=None).to_csv(out_dir / f'{species}_cpec_markers.csv', index=False)
data.write_h5ad(out_dir / f'{species}_cpec.h5ad')
atlas.obs['cell_subtype'] = atlas.obs.cell_type.astype(str)
atlas.obs.loc[data.obs_names, 'cell_subtype'] = data.obs.cell_type.astype(str)
atlas = atlas[~atlas.obs_names.isin(excluded)].copy()
atlas.obs['cell_subtype'] = atlas.obs.cell_subtype.astype('category')
atlas.write_h5ad(out_dir / f'{species}_final.h5ad')
