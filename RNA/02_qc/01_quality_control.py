#!/usr/bin/env python3
"""Apply species-specific QC and per-sample Scrublet."""

import os
from pathlib import Path
import re
import numpy as np
import pandas as pd
import scanpy as sc

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
samples = pd.read_csv(proj_dir / 'RNA/00_reference/samples.tsv', sep='\t', dtype=str)
out_dir = proj_dir / 'data/RNA/output/qc'
out_dir.mkdir(parents=True, exist_ok=True)
settings = {'human': (1000,8000,1000,30000,10,0.1,3), 'mouse': (500,4000,500,10000,10,0.15,2)}
counts = []
for species, limits in settings.items():
    sex_genes = set((proj_dir / f'RNA/00_reference/{species}_sex_genes.txt').read_text().splitlines())
    objects = []
    for row in samples.loc[(samples.include == '1') & (samples.species == species)].itertuples():
        data = sc.read_h5ad(proj_dir / f'data/RNA/output/load_data/{row.sample}.h5ad')
        n_input = data.n_obs
        if species == 'human':
            data = data[np.log10(data.obs.nuclear_counts + 1) > data.uns['nuclear_threshold']].copy()
        data.var['mt'] = data.var_names.str.upper().str.startswith('MT-')
        sc.pp.calculate_qc_metrics(data, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
        gmin,gmax,cmin,cmax,mt,ds,min_cells = limits
        keep = ((data.obs.n_genes_by_counts > gmin) & (data.obs.n_genes_by_counts < gmax) &
                (data.obs.total_counts > cmin) & (data.obs.total_counts < cmax) & (data.obs.pct_counts_mt < mt))
        data = data[keep].copy()
        sc.pp.scrublet(data, expected_doublet_rate=0.05, sim_doublet_ratio=2.0, random_state=0)
        data = data[data.obs.doublet_score < ds].copy()
        sc.pp.filter_genes(data, min_cells=min_cells)
        upper = data.var_names.str.upper()
        remove = (upper.str.startswith(('RPL','RPS')) | upper.str.match(r'^HB(?!P)') |
                  (upper == 'MALAT1') | data.var_names.isin(sex_genes))
        data = data[:, ~remove].copy()
        data.layers['counts'] = data.X.copy()
        objects.append(data)
        counts.append({'sample':row.sample, 'species':species, 'input':n_input, 'retained':data.n_obs})
    if objects:
        data = sc.concat(objects, join='inner', merge='same')
        data.write_h5ad(out_dir / f'{species}_counts.h5ad')
pd.DataFrame(counts).to_csv(out_dir / 'cell_counts.csv', index=False)
