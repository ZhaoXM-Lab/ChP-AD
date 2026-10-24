#!/usr/bin/env python3
"""Export cluster markers and annotate cell types."""

import os
from pathlib import Path
import pandas as pd
import scanpy as sc

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
species = os.environ.get('SPECIES', 'human')
out_dir = proj_dir / 'data/RNA/output/annotation'
out_dir.mkdir(parents=True, exist_ok=True)
data = sc.read_h5ad(proj_dir / f'data/RNA/output/umap/{species}_clustered.h5ad')
sc.tl.rank_genes_groups(data, groupby='cluster', method='wilcoxon', use_raw=True, pts=True)
sc.get.rank_genes_groups_df(data, group=None).to_csv(out_dir / f'{species}_cluster_markers.csv', index=False)
anno = pd.read_csv(proj_dir / 'RNA/00_reference/cluster_annotations.tsv', sep='\t', dtype=str)
anno = anno.loc[anno.species == species]
if anno.cluster.duplicated().any():
    raise ValueError('Each species/cluster annotation must be unique.')
mapping = anno.set_index('cluster').cell_type
data.obs['cell_type'] = data.obs.cluster.astype(str).map(mapping)
if data.obs.cell_type.isna().any():
    raise ValueError('Cluster annotations do not cover all observed clusters.')
data.obs['cell_type'] = data.obs.cell_type.astype('category')
data.write_h5ad(out_dir / f'{species}_annotated.h5ad')
