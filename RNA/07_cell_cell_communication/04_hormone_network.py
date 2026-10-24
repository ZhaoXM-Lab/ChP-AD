#!/usr/bin/env python3
"""Assemble the supplied sign-consistent hormone communication score."""

import os
from pathlib import Path
import numpy as np
import pandas as pd

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'data/RNA/output/hormone'
p = pd.read_csv(out / 'HPC_difference.csv').rename(columns={'cell_type':'sender'})
r = pd.read_csv(out / 'HRC_difference.csv').rename(columns={'cell_type':'receiver'})
edges = p.merge(r, on='Hormone', suffixes=('_producer','_receiver'))
edges = edges.loc[np.sign(edges.delta_producer) == np.sign(edges.delta_receiver)].copy()
edges['delta'] = (edges.delta_producer + edges.delta_receiver) / 2
edges = edges.loc[edges.delta != 0]
overlapping = (((edges.sender == 'CPEC') & edges.receiver.isin(['CPEC-1','CPEC-2'])) |
               ((edges.receiver == 'CPEC') & edges.sender.isin(['CPEC-1','CPEC-2'])))
edges = edges.loc[~overlapping].copy()
edges['mode'] = np.where(edges.sender == edges.receiver, 'Autocrine', 'Paracrine')
edges.to_csv(out / 'network_difference.csv', index=False)
