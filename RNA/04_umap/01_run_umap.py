#!/usr/bin/env python3
"""Construct the KNN graph, UMAP and Leiden clusters."""

import os
from pathlib import Path
import scanpy as sc

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
species = os.environ.get('SPECIES', 'human')
out_dir = proj_dir / 'data/RNA/output/umap'
out_dir.mkdir(parents=True, exist_ok=True)
data = sc.read_h5ad(proj_dir / f'data/RNA/output/integration/{species}_integrated.h5ad')
sc.pp.neighbors(data, use_rep='X_pca_harmony', n_neighbors=10, random_state=0)
sc.tl.umap(data, random_state=0)
sc.tl.leiden(data, resolution=0.5, key_added='cluster', random_state=0)
data.write_h5ad(out_dir / f'{species}_clustered.h5ad')
