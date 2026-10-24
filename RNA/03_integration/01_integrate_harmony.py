#!/usr/bin/env python3
"""Normalize, select 2,000 HVGs and integrate 30 PCs with Harmony."""

import os
from pathlib import Path
import scanpy as sc

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
species = os.environ.get('SPECIES', 'human')
out_dir = proj_dir / 'data/RNA/output/integration'
out_dir.mkdir(parents=True, exist_ok=True)
data = sc.read_h5ad(proj_dir / f'data/RNA/output/qc/{species}_counts.h5ad')
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
data.raw = data.copy()
sc.pp.highly_variable_genes(data, flavor='seurat', n_top_genes=2000)
data.var.loc[data.var.highly_variable].index.to_series().to_csv(out_dir / f'{species}_hvg.txt', index=False, header=False)
hvg = data[:, data.var.highly_variable].copy()
sc.pp.scale(hvg)
sc.pp.pca(hvg, n_comps=30, random_state=0)
sc.external.pp.harmony_integrate(hvg, key='sample', max_iter_harmony=30)
data.obsm['X_pca'] = hvg.obsm['X_pca']
data.obsm['X_pca_harmony'] = hvg.obsm['X_pca_harmony']
data.write_h5ad(out_dir / f'{species}_integrated.h5ad')
