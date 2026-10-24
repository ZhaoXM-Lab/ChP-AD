#!/usr/bin/env python3
"""Estimate hormone production and reception at major-lineage and subtype levels."""

import os
from pathlib import Path
import pandas as pd
import scanpy as sc
import hormone2cell as h2c

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'data/RNA/output/hormone'
out.mkdir(parents=True, exist_ok=True)
data = sc.read_h5ad(proj_dir / 'data/RNA/output/annotation/human_final.h5ad')
data = data[~data.obs.cell_type.isin(['Neuron','OPC','CD8+T'])].copy()
data.X = data.layers['counts'].copy()
os.chdir(out)
for condition in ['WT','AD']:
    combined = {'HPC': [], 'HRC': []}
    for level in ['cell_type','cell_subtype']:
        x = data[data.obs.disease == condition].copy()
        x.obs['suspension_type'] = 'nucleus'
        x.obs['Tissue'] = 'Brain'
        x.obs['hormone_cell_type'] = x.obs[level].astype(str)
        ave = h2c.compute_aveExp_by_category(x, sc_sn_col='suspension_type', celltype_col='hormone_cell_type', tissue_col='Tissue')
        ave = ave.loc[ave.suspension_type == 'nucleus']
        for role, definition, includes, excludes, adjust, threshold_file in [
            ('HPC',h2c.data.load_hormone_producing_file(),[f'hpc_include{i}' for i in range(1,9)],['hpc_exclude1','hpc_exclude2'],True,'max_expression_nucleus_pre.csv'),
            ('HRC',h2c.data.load_hormone_receptor_file(),['receptorgene1','receptorgene2'],[],False,'max_expression_nucleus_receptor2.csv')]:
            kwargs = dict(ave_all=ave, geneset_definition=definition, celltype_column='hormone_cell_type', tissue_col='Tissue',
                          adjustment=adjust, assay='nucleus', include_cols=includes, exclude_cols=excludes,
                          thresh_expr_low=0.01, thresh_pct=5, thresh_included_initial=0.1, thresh_excluded_initial=0.15,
                          use_precomputed='nucleus', max_expression_file=f'{condition}_{level}_{threshold_file}')
            if adjust:
                kwargs.update(specificity_threshold=0.8, coverage_threshold=0.5, thresh_included_adjusted=0.5, thresh_excluded_adjusted=0.15)
            result = h2c.hormone_strength(**kwargs).rename(columns={'hormone_cell_type':'cell_type'})
            if level == 'cell_subtype':
                result = result.loc[result.cell_type.isin(['CPEC-1','CPEC-2'])]
            combined[role].append(result)
    for role, frames in combined.items():
        pd.concat(frames, ignore_index=True).to_csv(out / f'{role}_{condition}.csv', index=False)
for role in ['HPC','HRC']:
    wt = pd.read_csv(out / f'{role}_WT.csv')
    ad = pd.read_csv(out / f'{role}_AD.csv')
    diff = wt[['Hormone','cell_type','Strength']].merge(ad[['Hormone','cell_type','Strength']], on=['Hormone','cell_type'], suffixes=('_WT','_AD'))
    diff['delta'] = diff.Strength_AD - diff.Strength_WT
    diff.to_csv(out / f'{role}_difference.csv', index=False)
