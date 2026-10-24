#!/usr/bin/env python3
"""Export raw counts, annotations and embeddings without dense conversion."""

import os
from pathlib import Path
import pandas as pd
import scanpy as sc
from scipy.io import mmwrite

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
species = os.environ.get('SPECIES', 'human')
data = sc.read_h5ad(proj_dir / f'data/RNA/output/annotation/{species}_final.h5ad')
out_dir = proj_dir / f'data/RNA/output/seurat_input/{species}'
out_dir.mkdir(parents=True, exist_ok=True)
mmwrite(str(out_dir / 'matrix.mtx'), data.layers['counts'].T, field='integer')
pd.Series(data.var_names).to_csv(out_dir / 'genes.tsv', sep='\t', header=False, index=False)
data.obs.to_csv(out_dir / 'cells.csv')
for key in ['X_pca_harmony','X_umap']:
    pd.DataFrame(data.obsm[key], index=data.obs_names).to_csv(out_dir / f'{key}.csv')
