#!/usr/bin/env python3
"""Load counts and identify human nuclei using the 30-transcript panel."""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.mixture import GaussianMixture
from skimage.filters import threshold_minimum

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out_dir = proj_dir / 'data/RNA/output/load_data'
out_dir.mkdir(parents=True, exist_ok=True)
samples = pd.read_csv(proj_dir / 'RNA/00_reference/samples.tsv', sep='\t', dtype=str)
nuclear_genes = ['MALAT1','NEAT1','FTX','FOXP1','RBMS3','ZBTB20','LRMDA','PBX1',
                 'ITPR2','AUTS2','TTC28','BNC2','EXOC4','RORA','PRKG1','ARID1B',
                 'PARD3B','GPHN','N4BP2L2','PKHD1L1','EXOC6B','FBXL7','MED13L','TBC1D5',
                 'IMMP2L','SYNE1','RERE','MBD5','EXT1','WWOX']

def nuclear_threshold(log_counts):
    gmm = GaussianMixture(n_components=2, means_init=[[1.0],[3.0]], tol=1, random_state=0)
    gmm.fit(log_counts.reshape(-1, 1))
    grid = np.linspace(0, log_counts.max(), 1000)
    density = np.exp(gmm.score_samples(grid[:, None]))
    sampled = np.repeat(grid, np.rint(density / density.sum() * 1e5).astype(int))
    return float(threshold_minimum(sampled))

if samples.empty:
    raise ValueError('RNA sample manifest is empty.')
for row in samples.loc[samples.include == '1'].itertuples():
    data = sc.read_10x_mtx(proj_dir / row.matrix_dir, var_names='gene_symbols', make_unique=False)
    if not data.var_names.is_unique:
        raise ValueError(f'{row.sample}: duplicate gene symbols in the input matrix.')
    data.obs['barcode'] = data.obs_names
    data.obs_names = data.obs_names + '-' + row.sample
    for key in ['sample','species','disease','source']:
        data.obs[key] = getattr(row, key)
    if row.species == 'human':
        panel = data.var_names.intersection(nuclear_genes)
        if len(panel) == 0:
            raise ValueError(f'{row.sample}: no nuclear-panel genes found.')
        data.obs['nuclear_counts'] = np.asarray(data[:, panel].X.sum(axis=1)).ravel()
        data.uns['nuclear_threshold'] = nuclear_threshold(np.log10(data.obs.nuclear_counts.to_numpy() + 1))
    data.write_h5ad(out_dir / f'{row.sample}.h5ad')
