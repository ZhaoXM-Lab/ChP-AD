#!/usr/bin/env python3
"""Summarize the same largest module across MCL inflation values."""

import os
from pathlib import Path
import pandas as pd

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
base = proj_dir / 'data/ATAC/output/grn'
selected = pd.read_csv(proj_dir / 'ATAC/00_reference/mcl_inflation.tsv', sep='\t')
stats = []
for folder in sorted((base / 'mcl').iterdir()):
    for file in sorted(folder.glob('Clusters_I*.txt')):
        groups = [line.split() for line in file.read_text().splitlines() if line.strip()]
        largest = max(groups, key=len)
        stats.append({'cell_type':folder.name,'inflation':float(file.stem.split('_I')[1]),'cluster_count':len(groups),
                      'largest_size':len(largest),'largest_TFs':sum(x.startswith('TF_') for x in largest),
                      'largest_genes':sum(x.startswith('G_') for x in largest)})
pd.DataFrame(stats).to_csv(base / 'mcl_parameter_statistics.csv', index=False)
all_nodes = []
module_stats = []
for row in selected.itertuples():
    file = base / 'mcl' / row.cell_type / f'Clusters_I{row.inflation:.1f}.txt'
    edges = pd.read_csv(base / f'{row.cell_type}_grn.csv')
    for i, line in enumerate(file.read_text().splitlines(), 1):
        nodes = line.split()
        targets = {x[2:] for x in nodes if x.startswith('G_')}
        tfs = {x[3:] for x in nodes if x.startswith('TF_')}
        deg = set(edges.loc[edges.direction != 'NS','gene'])
        module_stats.append({'cell_type':row.cell_type,'module':i,'n_TFs':len(tfs),'n_targets':len(targets),
                             'n_DEGs':len(targets & deg),'DEG_ratio':len(targets & deg)/len(targets) if targets else 0})
        for node in nodes:
            all_nodes.append({'cell_type':row.cell_type,'module':i,'node':node,'name':node[3:] if node.startswith('TF_') else node[2:],
                              'type':'TF' if node.startswith('TF_') else 'Gene'})
pd.DataFrame(all_nodes).to_csv(base / 'module_nodes.csv', index=False)
pd.DataFrame(module_stats).sort_values(['cell_type','DEG_ratio'], ascending=[True,False]).to_csv(base / 'module_statistics.csv', index=False)
