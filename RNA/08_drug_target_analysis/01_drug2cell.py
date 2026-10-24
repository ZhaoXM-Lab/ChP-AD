#!/usr/bin/env python3
"""Score drug targets and intersect them with cell-type-specific AD DEGs."""

import os
from pathlib import Path
import pandas as pd
import scanpy as sc
import drug2cell as d2c

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'data/RNA/output/drug2cell'
out.mkdir(parents=True, exist_ok=True)
data = sc.read_h5ad(proj_dir / 'data/RNA/output/annotation/human_final.h5ad')
data = data[~data.obs.cell_type.isin(['Neuron','OPC','CD8+T'])].copy()
data.X = data.layers['counts'].copy()
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
deg = pd.read_csv(proj_dir / 'data/RNA/output/differential_expression/consensus_degs.csv')
outputs = []
for level in ['cell_type','cell_subtype']:
    x = data.copy()
    x.obs['drug_group'] = x.obs[level].astype('category')
    sc.tl.rank_genes_groups(x, groupby='drug_group', method='wilcoxon', use_raw=False)
    markers = sc.get.rank_genes_groups_df(x, group=None)
    markers = markers.loc[(markers.pvals_adj < 0.05) & (markers.logfoldchanges > 1)]
    d2c.score(x, nested=True, method='mean')
    scores = x.uns['drug2cell']
    sc.tl.rank_genes_groups(scores, groupby='drug_group', method='wilcoxon', pts=True)
    rank = sc.get.rank_genes_groups_df(scores, group=None)
    rank = rank.loc[(rank.pvals_adj < 0.05) & (rank.logfoldchanges > 1)].copy()
    targets = d2c.util.prepare_targets(x, targets=None, reconstruct=False)
    for ct in rank.group.unique():
        if level == 'cell_subtype' and ct not in ['CPEC-1','CPEC-2']:
            continue
        disease_genes = set(deg.loc[(deg.cell_type == ct) & deg.specific, 'gene'])
        marker_genes = set(markers.loc[markers.group == ct, 'names'])
        eligible = disease_genes & marker_genes
        for row in rank.loc[rank.group == ct].itertuples():
            for gene in sorted(set(targets.get(row.names, [])) & eligible):
                outputs.append({'drug':row.names,'gene':gene,'cell_type':ct,'score':row.scores,'p_adj':row.pvals_adj})
    rank.to_csv(out / f'{level}_drug_scores.csv', index=False)
pd.DataFrame(outputs, columns=['drug','gene','cell_type','score','p_adj']).to_csv(out / 'drug_gene_interactions.csv', index=False)
