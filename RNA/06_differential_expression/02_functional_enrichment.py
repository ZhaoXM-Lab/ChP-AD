#!/usr/bin/env python3
"""Run GO-BP Enrichr and CPEC preranked GSEA."""

import os
from pathlib import Path
import pandas as pd
import gseapy as gp

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'data/RNA/output/enrichment'
out.mkdir(parents=True, exist_ok=True)
library = 'GO_Biological_Process_2025'
deg = pd.read_csv(proj_dir / 'data/RNA/output/differential_expression/consensus_degs.csv')
for (ct, direction), group in deg.loc[deg.specific & (deg.cell_type != 'Pericyte')].groupby(['cell_type','direction']):
    gp.enrichr(gene_list=group.gene.tolist(), gene_sets=library, organism='human',
               outdir=str(out / f'{ct}_{direction}'), cutoff=0.05)
markers = pd.read_csv(proj_dir / 'data/RNA/output/annotation/human_cpec_markers.csv')
for ct, group in markers.groupby('group'):
    rank = group[['names','logfoldchanges']].dropna().drop_duplicates('names').sort_values('logfoldchanges', ascending=False)
    gp.prerank(rnk=rank, gene_sets=library, min_size=15, max_size=500,
               permutation_num=1000, seed=123, threads=4, outdir=str(out / f'{ct}_gsea'))
    selected = group.loc[(group.pvals_adj < 0.05) & (group.logfoldchanges > 0.3), 'names']
    gp.enrichr(gene_list=selected.tolist(), gene_sets=library, organism='human', outdir=str(out / f'{ct}_markers'), cutoff=0.05)
