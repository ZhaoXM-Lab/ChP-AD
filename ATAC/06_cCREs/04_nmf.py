#!/usr/bin/env python3
"""Factor chromatin accessibility into eight modules."""

import os
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from scipy.special import xlogy

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
out = proj_dir / 'data/ATAC/output/ccres/nmf'
data = pd.read_csv(out / 'accessibility.tsv', sep='\t', index_col=0)
if np.any(data.to_numpy() < 0) or not np.isfinite(data.to_numpy()).all():
    raise ValueError('NMF requires finite, nonnegative accessibility.')
model = NMF(n_components=8, init='nndsvda', max_iter=1000, random_state=1)
w = model.fit_transform(data)
h = model.components_
modules = [f'M{i+1}' for i in range(8)]
pd.DataFrame(w, index=data.index, columns=modules).to_csv(out / 'W.tsv', sep='\t')
pd.DataFrame(h, index=modules, columns=data.columns).to_csv(out / 'H.tsv', sep='\t')
norm = np.divide(w, w.sum(axis=1, keepdims=True), out=np.zeros_like(w), where=w.sum(axis=1, keepdims=True) != 0)
specificity = 1 + xlogy(norm, norm).sum(axis=1) / np.log(8)
specificity[w.sum(axis=1) == 0] = 0
median = np.median(specificity)
mad = np.median(np.abs(specificity - median))
keep = (specificity > median + mad) & (norm.max(axis=1) > np.median(norm.max(axis=1)))
pd.DataFrame({'peak':data.index[keep], 'module':np.array(modules)[w[keep].argmax(axis=1)],
              'specificity':specificity[keep]}).to_csv(out / 'module_peaks.csv', index=False)
activity = pd.DataFrame(h / h.max(), index=modules, columns=data.columns)
activity.to_csv(out / 'normalized_activity.tsv', sep='\t')
active = activity.stack().rename('activity').reset_index()
active.columns = ['module','cell_type_diagnosis','activity']
active.loc[active.activity > np.quantile(activity, 0.95)].to_csv(out / 'active_modules.csv', index=False)
trends = []
for k in range(2, min(20, min(data.shape)) + 1):
    fit = NMF(n_components=k, init='nndsvda', max_iter=1000, random_state=1)
    fit.fit_transform(data)
    nh = fit.components_ / np.maximum(fit.components_.sum(axis=0, keepdims=True), 1e-12)
    sparse = (np.sqrt((nh**2).sum(axis=0)) - 1/np.sqrt(k)) / (1 - 1/np.sqrt(k))
    trends.append({'k':k,'mean_sparsity':sparse.mean(),'mean_entropy':(-xlogy(nh,nh).sum(axis=0)/np.log(2)).mean()})
pd.DataFrame(trends).to_csv(out / 'rank_diagnostics.csv', index=False)
