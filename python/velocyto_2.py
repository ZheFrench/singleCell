import numpy as np
import pandas as pd
import anndata
from scanpy import read
import argparse
import os
import scvelo as scv
import igraph

# First 
# loom_combine.R
# export HDF5_USE_FILE_LOCKING='FALSE'
# ssh -X
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", action="store", help="file h5ad", required=True, type=str, dest='h5ad')
parameters = parser.parse_args()

print(parameters.h5ad)
file = os.path.basename(parameters.h5ad)
print(file)
path = os.path.abspath(parameters.h5ad)
print(path)

sample_one = scv.read(parameters.h5ad)

###############################################################
# scvelo
###############################################################
print ("scvelo")
scv.pl.proportions(sample_one,save = 'proportion.png',show=False)

scv.pp.filter_and_normalize(sample_one, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(sample_one, n_pcs=30, n_neighbors=30)
scv.tl.velocity(sample_one)
scv.tl.velocity_graph(sample_one)

# Project the velocities 
print(sample_one.uns)
print(sample_one.obsm)

#odict_keys(['pca', 'neighbors', 'velocity_params', 'velocity_graph', 'velocity_graph_neg'])

#KeysView(AxisArrays with keys: X_umap, X_pca)

scv.pl.umap(sample_one,  frameon=False, legend_loc='on data', color="seurat_clusters",title='', save=path+'.celltypes2.png',show=False)

scv.pl.velocity_embedding_stream(sample_one, basis="umap", color="seurat_clusters",save = path+'.stream2.png',show=False)
scv.pl.velocity_embedding_stream(sample_one, basis="tsne", color="seurat_clusters",save = path+'.streamTsne2.png',show=False)

# Scatter
scv.pl.velocity_embedding(sample_one, basis="umap", color="seurat_clusters", arrow_length=3, arrow_size=2, dpi=120,save = path+'.scatter2.png',show=False)

scv.tl.recover_dynamics(sample_one)
scv.tl.latent_time(sample_one)
scv.pl.scatter(sample_one, color="latent_time", color_map="gnuplot",save = path+'.latentTime2.png',show=False)

scv.tl.velocity_confidence(sample_one)#s_genes g2m_genes
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(sample_one, c=keys, cmap='coolwarm', perc=[5, 95],save = path+'.speed_and_confidence2.png',show=False)

#WARNING: No root cells detected. Consider specifying root cells to improve latent time prediction.
#computing latent time using root_cells as prior
#No root cells detected

#top_genes = sample_one.var["fit_likelihood"].sort_values(ascending=False).index[:300]
#scv.pl.heatmap(sample_one, var_names=top_genes, sortby="latent_time", col_color="seurat_clusters", n_convolve=100,save = file+'heatmap.png',show=False)

scv.pl.velocity(sample_one, var_names=['MKI67','ANOS1','CRYAB','AXL'],save = path+'.genes2.png', nrows=5, show=False)

scv.pl.scatter(sample_one, 'HES4', color=['seurat_clusters', 'velocity'],save = path+'.genes_selected2.png',   add_outline=True, show=False)

scv.tl.rank_velocity_genes(sample_one, groupby='seurat_clusters', min_corr=.3)

df = scv.DataFrame(sample_one.uns['rank_velocity_genes']['names'])
#df.head()

# this is needed due to a current bug - bugfix is coming soon.
sample_one.uns['neighbors']['distances'] = sample_one.obsp['distances']
sample_one.uns['neighbors']['connectivities'] = sample_one.obsp['connectivities']

scv.tl.paga(sample_one, groups='seurat_clusters')
#df = scv.get_df(sample_one, 'paga/transitions_confidence', precision=2).T
#df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(sample_one, basis='umap', color="seurat_clusters", size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5,save = path+'.paga2.png',   add_outline=True, show=False)

#tempsample_one = sample_one.raw.to_sample_one()
#tempsample_one.var_names = sample_one.var_names
#sample_one.raw = tempsample_one
#scv.tl.score_genes_cell_cycle(sample_one)
#scv.pl.scatter(sample_one.raw, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95],save = 'cycling.png',show=False)

#######################################################################################
