import sys
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import cell2location
import scvi
import os
from matplotlib import rcParams
rcParams['pdf.fonttype'] = 42
import seaborn as sns
import time
import gc
gc.collect()


# load data
cellcount = "/home/myi/data/scRNA-seq/sc_counts_SHK166_RA_Knee.mtx"
sp_data = "/home/myi/data/spatial-trancriptomics/st_counts_SHK166_RA_Knee.mtx"
celltype = "/home/myi/data/scRNA-seq/ctypes7_counts_SHK166_RA_Knee.csv"

start_time = time.time()
adata_ref = sc.read_mtx(cellcount).transpose()
df_celltype = pd.read_csv(celltype, header=0, index_col=0)
df_celltype.index = adata_ref.obs.index
adata_ref.obs['Subset'] = df_celltype

shuffled_obs_indices = np.random.permutation(adata_ref.n_obs)
adata_ref = adata_ref[shuffled_obs_indices, :]

# Get 'Method'
# adata_ref.obs['Method'] = '3GEX'

# Get 'Sample'
# adata_ref.obs['Sample'] = adata_ref.obs_names
# adata_ref.obs['Sample'] = adata_ref.obs['Sample'].apply(lambda x: x[0:4])
N_cells_per_location = 7
# QC
from cell2location.utils.filtering import filter_genes
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
# 5, 0.03, 1.12

# In our case, a few genes are cut
adata_ref = adata_ref[:, selected].copy()

# setup scvi
cell2location.models.RegressionModel.setup_anndata(adata=adata_ref,
                        # 10X reaction / sample / batch
                        # batch_key='Sample',
                        # cell type, covariate used for constructing signatures
                        labels_key='Subset',
                        # multiplicative technical effects (platform, 3' vs 5', donor effect)
                        # categorical_covariate_keys=['Method']
                       )

from cell2location.models import RegressionModel
mod = RegressionModel(adata_ref)
mod.view_anndata_setup(adata_ref)

# Use all data for training (validation not implemented yet, train_size=1)
mod.train(max_epochs=250, batch_size=2500, train_size=1, lr=0.002, use_gpu=True)

# plot ELBO loss history during training, removing first 20 epochs from the plot
mod.plot_history(20)
plt.savefig(f"./history.png",
                bbox_inches='tight')
plt.close()

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_ref = mod.export_posterior(
    adata_ref, sample_kwargs={'num_samples': 1000, 'batch_size': 2500, 'use_gpu': True}
)

# Results saving folder
ref_run_name = f'./results/analysis/reference_signatures'

# Save model
mod.save(f"{ref_run_name}", overwrite=True)

adata_file = f"{ref_run_name}/sc.h5ad"
adata_ref.write(adata_file)

# First plot: Reconstruction accuracy to assess if there are any issues with inference
# Second plot: When this plot is very different from a diagonal plot (e.g. very low values on Y-axis, density everywhere) it indicates problems with signature estimation.
mod.plot_QC()
plt.savefig(f"./reconstruction_accuracy_histogram.png",
                bbox_inches='tight')
plt.close()

# Export estimated expression in each cluster
if 'means_per_cluster_mu_fg' in adata_ref.varm.keys():
    inf_aver = adata_ref.varm['means_per_cluster_mu_fg'][[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
else:
    inf_aver = adata_ref.var[[f'means_per_cluster_mu_fg_{i}'
                                    for i in adata_ref.uns['mod']['factor_names']]].copy()
inf_aver.columns = adata_ref.uns['mod']['factor_names']
inf_aver.iloc[0:5, 0:5]

# Note that the column names are cell types, the row names are gene names, in the original code
# the authors use ENSEMBL as names instead of raw gene names

# ST
adata_vis = sc.read_mtx(sp_data).transpose() # cells by genes

adata_vis.obs['sample'] = 'SHK' # Since it is manually generated

# find shared genes and subset both anndata and reference signatures
intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
adata_vis = adata_vis[:, intersect].copy()
inf_aver = inf_aver.loc[intersect, :].copy()

# prepare anndata for cell2location model
cell2location.models.Cell2location.setup_anndata(adata=adata_vis, batch_key="sample")

# create and train the model
mod = cell2location.models.Cell2location(
    adata_vis, cell_state_df=inf_aver,
    # the expected average cell abundance: tissue-dependent
    # hyper-prior which can be estimated from paired histology:
    N_cells_per_location=N_cells_per_location,
    # hyperparameter controlling normalisation of
    # within-experiment variation in RNA detection:
    detection_alpha=2
)
mod.view_anndata_setup()

mod.train(max_epochs=10000,
          # train using full data (batch_size=None)
          batch_size=None,
          # use all data points in training because
          # we need to estimate cell abundance at all locations
          train_size=1,
          use_gpu=True,
         )

# plot ELBO loss history during training, removing first 100 epochs from the plot
mod.plot_history(1000)
plt.legend(labels=['full data training']);

# In this section, we export the estimated cell abundance (summary of the posterior distribution).
adata_vis = mod.export_posterior(
    adata_vis, sample_kwargs={'num_samples': 1000, 'batch_size': mod.adata.n_obs, 'use_gpu': True}
)

mod.plot_QC()
plt.savefig(f"./posterior.png",
                bbox_inches='tight')
plt.close()

# Save model
run_name = "cell2loc"
mod.save(f"{run_name}", overwrite=True)

# mod = cell2location.models.Cell2location.load(f"{run_name}", adata_vis)

# Save anndata object with results
adata_file = f"{run_name}/sp.h5ad"
adata_vis.write(adata_file)
adata_file
