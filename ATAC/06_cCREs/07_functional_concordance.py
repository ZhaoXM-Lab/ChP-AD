#!/usr/bin/env python3
"""Compare the top 100 GO-BP terms from GREAT and RNA Enrichr."""

import os
from pathlib import Path
import pandas as pd

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
base = proj_dir / 'data/ATAC/output/ccres'
rows = []
markers = pd.read_csv(base / 'specific_ccres.csv')
for ct in sorted(set(markers.cell_type) - {'Pericyte','Neuron','OPC','CD8+T'}):
    for disease,direction in [('AD','Up'),('WT','Down')]:
        atac = pd.read_csv(base / f'great/{ct}_{disease}.csv')
        rna = pd.read_csv(proj_dir / f'data/RNA/output/enrichment/{ct}_{direction}/GO_Biological_Process_2025.human.enrichr.reports.txt',sep='\t')
        a = set(atac.sort_values('p_adjust_hyper').head(100).id.dropna())
        b = set(rna.sort_values('Adjusted P-value').head(100).Term.str.extract(r'(GO:\d+)',expand=False).dropna())
        common = a & b
        rows.append({'cell_type':ct,'disease':disease,'n_GREAT':len(a),'n_Enrichr':len(b),'n_overlap':len(common),
                     'overlap_ratio':0.5*(len(common)/len(a)+len(common)/len(b)) if min(len(a),len(b)) else float('nan')})
pd.DataFrame(rows).to_csv(base / 'functional_concordance.csv',index=False)
