#!/usr/bin/env python3
"""Fit CellOracle and rank virtual TF knockouts using pseudotime perturbation scores."""

import os
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import celloracle as co
from celloracle.applications import Oracle_development_module, Gradient_calculator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

proj_dir = Path(os.environ['PROJ_DIR']).resolve()
ct = os.environ.get('CELL_TYPE','Macrophage')
out = proj_dir / 'data/ATAC/output/celloracle' / ct
data = sc.read_h5ad(proj_dir / 'data/RNA/output/annotation/human_final.h5ad')
column = 'cell_subtype' if ct in ['CPEC-1','CPEC-2'] else 'cell_type'
data = data[data.obs[column] == ct].copy()
trajectory = pd.read_csv(out / 'trajectory.csv').set_index('cell', verify_integrity=True)
if set(data.obs_names) != set(trajectory.index):
    raise ValueError('RNA and Monocle3 cell identities differ.')
data = data[trajectory.index].copy()
if not np.isfinite(trajectory.pseudotime).all():
    raise ValueError('Monocle3 returned nonfinite pseudotime.')
hvg = (out / 'hvg.txt').read_text().splitlines()
data = data[:, data.var_names.isin(hvg)].copy()
data.X = data.layers['counts'].copy()
data.layers['raw_count'] = data.X.copy()
sc.pp.normalize_total(data, target_sum=1e4)
sc.pp.log1p(data)
data.obs['oracle_cluster'] = ct
data.obs['monocle3_pseudotime'] = trajectory.pseudotime
data.obsm['X_umap'] = trajectory[['UMAP1','UMAP2']].to_numpy()
grn = pd.read_csv(proj_dir / f'data/ATAC/output/grn/{ct}_grn.csv')
grn = grn.loc[grn.TF.isin(data.var_names) & grn.gene.isin(data.var_names)]
pd.Series(sorted(grn.TF.unique()), name='TF').to_csv(out / 'input_HVG_TFs.csv', index=False)
base = (pd.crosstab(grn.gene, grn.TF) > 0).astype(int)
base['gene_short_name'] = base.index
base['peak_id'] = 'integrated_grn'
oracle = co.Oracle()
oracle.import_anndata_as_normalized_count(adata=data, cluster_column_name='oracle_cluster', embedding_name='X_umap')
oracle.import_TF_data(TF_info_matrix=base)
oracle.perform_PCA()
k = max(5, int(0.025 * data.n_obs))
oracle.knn_imputation(n_pca_dims=30, k=k, balanced=True, b_sight=k*8, b_maxl=k*4, n_jobs=4)
links = oracle.get_links(cluster_name_for_GRN_unit='oracle_cluster', alpha=10, verbose_level=1)
links.filter_links(p=0.01, weight='coef_abs', threshold_number=2000)
links.links_dict[ct].to_csv(out / 'fitted_links.csv', index=False)
oracle.get_cluster_specific_TFdict_from_Links(links_object=links)
oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)
oracle.to_hdf5(str(out / f'{ct}.celloracle.oracle'))
oracle.adata.obs['Pseudotime'] = oracle.adata.obs.monocle3_pseudotime
gradient = Gradient_calculator(oracle_object=oracle)
gradient.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
gradient.calculate_mass_filter(min_mass=0.01)
gradient.transfer_data_into_grid(args={'method':'polynomial','n_poly':3})
gradient.calculate_gradient()
dev = Oracle_development_module()
dev.load_differentiation_reference_data(gradient_object=gradient)
results = []
example_tf = {'CPEC':'TCF4','Macrophage':'AHR'}.get(ct)
if example_tf and example_tf not in set(links.links_dict[ct].source):
    raise ValueError(f'{example_tf} is absent from the fitted {ct} network.')
for tf in sorted(links.links_dict[ct].source.unique()):
    np.random.seed(0)
    oracle._clear_simulation_results()
    oracle.simulate_shift(perturb_condition={tf:0.0}, n_propagation=3, ignore_warning=True)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)
    dev.load_perturb_simulation_data(oracle_object=oracle)
    dev.calculate_inner_product()
    dev.calculate_digitized_ip(n_bins=10)
    score = np.asarray(dev.inner_product)
    results.append({'TF':tf,'mean_PS':score.mean(),'sum_PS':score.sum(),'sum_neg_PS':score[score < 0].sum()})
    if tf == {'CPEC':'TCF4','Macrophage':'AHR'}.get(ct):
        oracle.calculate_p_mass(smooth=0.8, n_grid=40, n_neighbors=200)
        oracle.calculate_mass_filter(min_mass=0.01)
        fig, ax = plt.subplots(figsize=(6,6))
        coords = np.asarray(oracle.embedding)
        for d, color in [('AD','#cc0000'),('WT','#00167a')]:
            mask = (oracle.adata.obs.disease == d).to_numpy()
            ax.scatter(coords[mask,0],coords[mask,1],c=color,s=3,alpha=0.25,rasterized=True)
        oracle.plot_simulation_flow_on_grid(scale=25, ax=ax, show_background=False)
        ax.set_axis_off()
        fig.savefig(out / f'{tf}_KO.pdf', bbox_inches='tight')
        plt.close(fig)
pd.DataFrame(results).sort_values('sum_neg_PS').to_csv(out / 'TF_perturb_score_by_inhibition.csv', index=False)
