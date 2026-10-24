#!/usr/bin/env python3
"""Collect enrichment and one-sided coefficient P values without filling missing results."""

import os
from pathlib import Path
import pandas as pd
from math import erfc, sqrt

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
base = proj_dir / 'data/ATAC/output/ldsc'
rows = []
for file in sorted((base / 'results').glob('*/*.results')):
    table = pd.read_csv(file, sep=r'\s+')
    selected = table.loc[table.Category == 'ANNOTL2_0']
    if len(selected) != 1:
        raise ValueError(f'{file.name}: expected the ANNOTL2_0 custom annotation.')
    row = selected.iloc[0]
    z = row['Coefficient'] / row['Coefficient_std_error']
    rows.append({'cell_type':file.parent.name,'trait':file.stem,'enrichment':row['Enrichment'],
                 'coefficient':row['Coefficient'],'coefficient_se':row['Coefficient_std_error'],
                 'p_one_sided':0.5 * erfc(z / sqrt(2))})
if not rows:
    raise ValueError('No LDSC regression results found.')
pd.DataFrame(rows).to_csv(base / 'enrichment_summary.csv', index=False)
