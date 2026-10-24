#!/usr/bin/env python3
"""Plot enrichment z-scores, coefficient significance and correlations with AD."""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'figures/ATAC'
out.mkdir(parents=True,exist_ok=True)
data = pd.read_csv(proj_dir / 'data/ATAC/output/ldsc/enrichment_summary.csv')
order = pd.read_csv(proj_dir / 'ATAC/00_reference/gwas.tsv',sep='\t').trait.tolist()
matrix = data.pivot(index='cell_type',columns='trait',values='enrichment').reindex(columns=order)
z = (matrix-matrix.mean())/matrix.std(ddof=0)
p = data.pivot(index='cell_type',columns='trait',values='p_one_sided').reindex(index=z.index,columns=z.columns)
stars = p.map(lambda x: '***' if x < 0.01 else '**' if x < 0.05 else '*' if x < 0.1 else '')
fig,ax = plt.subplots(figsize=(8,5))
sns.heatmap(z,cmap='RdBu_r',center=0,annot=stars,fmt='',ax=ax)
fig.savefig(out / 'Fig6a_LDSC.pdf',bbox_inches='tight')
plt.close(fig)
correlations = matrix.corr(method='pearson')['AD'].drop('AD').sort_values()
ax = correlations.plot.barh(figsize=(4,4))
ax.set_xlabel('Pearson r with AD enrichment')
ax.figure.savefig(out / 'Fig6b_AD_correlations.pdf',bbox_inches='tight')
plt.close(ax.figure)
