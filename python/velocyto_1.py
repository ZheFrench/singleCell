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
parser.add_argument("-f", "--file", action="store", help="file h5ad", required=True, type=str, dest='loom'),
parser.add_argument("-c", "--conditions", action="store", help="file h5ad", required=True, type=str, dest='conditions'),
parser.add_argument("-d", "--dir", action="store", help="file h5ad", required=True, type=str, dest='dir')

parameters = parser.parse_args()

print(parameters.loom)
file = os.path.basename(parameters.loom)
path = os.path.abspath(parameters.loom)

print(file)
print(path)

dir = parameters.dir 
conditions = parameters.conditions
###############################################################
# Metsample_one
###############################################################

print("umap_cord")
umap_cord     = pd.read_csv(dir+conditions+"_cell_embeddings.csv")

umap_cord.rename(columns={'Unnamed: 0':'CellID'}, inplace=True )

umap_cord['CellID'] = umap_cord['CellID'].str[::-1]
umap_cord['CellID'] = umap_cord['CellID'].str.replace('_',':',1)
umap_cord['CellID'] = umap_cord['CellID'].str[::-1]
umap_cord['CellID'] = umap_cord['CellID'].str.replace('-1','x')
print(umap_cord.head())

print("sample_obs")
sample_obs    = pd.read_csv(dir+conditions+"_cellID_obs.csv")
sample_obs.rename(columns={'Cells(seurat.Object)':'CellID'}, inplace=True )

sample_obs['CellID'] = sample_obs['CellID'].str[::-1]
sample_obs['CellID'] = sample_obs['CellID'].str.replace('_',':',1)
sample_obs['CellID'] = sample_obs['CellID'].str[::-1]
sample_obs['CellID'] = sample_obs['CellID'].str.replace('-1','x')

print(sample_obs.head())

print("cell_clusters")
cell_clusters = pd.read_csv(dir+conditions+"_clusters.csv")

cell_clusters['CellID'] = cell_clusters['CellID'].str[::-1]
cell_clusters['CellID'] = cell_clusters['CellID'].str.replace('_',':',1)
cell_clusters['CellID'] = cell_clusters['CellID'].str[::-1]
cell_clusters['CellID'] = cell_clusters['CellID'].str.replace('-1','x')

print(cell_clusters.head())


print("ident")
ident = pd.read_csv(dir+conditions+"_ident.csv")

ident['CellID'] = ident['CellID'].str[::-1]
ident['CellID'] = ident['CellID'].str.replace('_',':',1)
ident['CellID'] = ident['CellID'].str[::-1]
ident['CellID'] = ident['CellID'].str.replace('-1','x')

help_dict = {
    'CTL_Vertes': '#7DE9BF',
    'CTL_Rouges': '#FF876E',
    '4006_rouge': '#922008',
    '4006_verte': '#129360',
     'OSI_TIPI_Rouges': '#FF1300',
    'OSI_TIPI_Vertes': '#82FC5B', 
    'Human_CTL_GRCh38' :'#2E4053',
    'Human_OSI_TIPI' :'#D35400',
    'Human_OSI' :'#C0392B',
     'Human_TIPI':'#8E44AD'

}
ident = ident.replace({"orig.ident": help_dict})

print(ident.head())




###############################################################
# Loop
###############################################################

print("Loop")
# Loom files were created from the directory ouput of 10X
# They were combined with loom_combine.py
# They need to be filtered as we did previously for each sample. ( we used stuffs exported from plotNestedCluster.R)
sample_one = anndata.read_loom(parameters.loom)
#sample_one.var_names_make_unique()
sample_one.var_names_make_unique()

sample_one = sample_one[np.isin(sample_one.obs.index,sample_obs["CellID"])]
sample_one_index = pd.DataFrame(sample_one.obs.index)

print("cell_clusters_ordered")
#"Now if we merge our index dataframe with our UMAP, the order will match our anndata object."
cell_clusters_ordered = sample_one_index.merge(cell_clusters, on = "CellID")
#"Since we're certain the orders are the same, we can remove the first column of the data frame and add the UMAP coordinates to our anndata object."
cell_clusters_ordered = cell_clusters_ordered.iloc[:,1:]
sample_one.uns['cluster'] = cell_clusters_ordered.values
print(cell_clusters_ordered)

print("ident_ordered")
ident_ordered = sample_one_index.merge(ident, on = "CellID")
ident_ordered = ident_ordered.iloc[:,1:]
sample_one.uns['ident'] = ident_ordered.values
print(ident_ordered)

print("umap_ordered")
umap_ordered = sample_one_index.merge(umap_cord, on = "CellID")
umap_ordered = umap_ordered.iloc[:,1:]
sample_one.obsm['X_umap'] = umap_ordered.values
print(umap_ordered)

###############################################################
# scvelo
###############################################################
print ("scvelo")
scv.pl.proportions(sample_one,save = path+'.proportion.png',show=False)

scv.pp.filter_and_normalize(sample_one, min_shared_counts = 20, n_top_genes = 2000) #, min_shared_counts = 20, n_top_genes = 2000
scv.pp.moments(sample_one, n_pcs = 30, n_neighbors = 30) #
scv.tl.velocity(sample_one)
scv.tl.velocity_graph(sample_one)

# Project the velocities 
print("uns")
print(sample_one.uns.keys())
print("obsm")
print(sample_one.obsm.keys())


#seurat_clusters
#scv.pl.umap(sample_one,  frameon = False, color  = sample_one.uns['cluster'], title='',legend_loc='right', save=path+'.celltypes.png',show=False)

#scv.pl.velocity_embedding_stream(sample_one, basis="umap", color  = sample_one.uns['cluster'],legend_loc='right',save = path+'.cluster.stream.png',show=False)
#scv.pl.velocity_embedding_stream(sample_one, basis="umap", color  = sample_one.uns['ident'],legend_loc='right',save = path+'.ident.stream.png',show=False)

#scv.pl.velocity_embedding_grid(sample_one, basis="umap", color  = sample_one.uns['cluster'],legend_loc='right',arrow_length=6, arrow_size=4,save = path+'.grid.png',show=False)

# Scatter
#scv.pl.velocity_embedding(sample_one, basis="umap", color  = sample_one.uns['cluster'], legend_loc='right',arrow_length=6, arrow_size=4, dpi=120,save = path+'.cluster.scatter.png',show=False)
#scv.pl.velocity_embedding(sample_one, basis="umap", color  = sample_one.uns['ident'], legend_loc='right',arrow_length=6, arrow_size=4, dpi=120,save = path+'.ident.scatter.png',show=False)

##################################
#Testing top-likelihood genes
##################################

scv.tl.recover_dynamics(sample_one)
scv.tl.latent_time(sample_one)
scv.pl.scatter(sample_one, color="latent_time", color_map="gnuplot",save = path+'.latentTime.png',show=False)

scv.tl.velocity_confidence(sample_one)#s_genes g2m_genes
keys = 'velocity_length', 'velocity_confidence'
#scv.pl.scatter(sample_one, c=keys, cmap='coolwarm', perc=[5, 95],save = path+'.speed_and_confidence.png',show=False)
scv.tl.rank_velocity_genes(sample_one, groupby=sample_one.uns['cluster'], min_corr=.3)
df = scv.DataFrame(sample_one.uns['rank_velocity_genes']['names'])
df.write( path+'.rank_velocity_genes.csv')
#WARNING: No root cells detected. Consider specifying root cells to improve latent time prediction.
#computing latent time using root_cells as prior
#No root cells detected

top_genes = sample_one.var["fit_likelihood"].sort_values(ascending=False).index[:100]
scv.tl.differential_kinetic_test(sample_one, var_names=top_genes, groupby=sample_one.uns['cluster'])

scv.pl.heatmap(sample_one, var_names=top_genes, sortby="latent_time", col_color=sample_one.uns['cluster'], n_convolve=100,save = file+'heatmap.png',show=False)
#adata.write('data/pancreas.h5ad')
df = scv.DataFrame(top_genes)
df.write( path+'.top_genes.csv')

exit(0)

scv.pl.velocity(sample_one, var_names=['MKI67','ANOS1','CRYAB','AXL'],save = path+'.genes.png', nrows=5, show=False)

scv.pl.scatter(sample_one, 'HES4', color=[sample_one.uns['cluster'], 'velocity'],save = path+'.genes_selected.png',   add_outline=True, show=False)

scv.tl.rank_velocity_genes(sample_one, groupby=sample_one.uns['cluster'], min_corr=.3)

df = scv.DataFrame(sample_one.uns['rank_velocity_genes']['names'])
#df.head()

# this is needed due to a current bug - bugfix is coming soon.
sample_one.uns['neighbors']['distances'] = sample_one.obsp['distances']
sample_one.uns['neighbors']['connectivities'] = sample_one.obsp['connectivities']

scv.tl.paga(sample_one, groups=sample_one.uns['cluster'])
#df = scv.get_df(sample_one, 'paga/transitions_confidence', precision=2).T
#df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(sample_one, basis='umap', color=sample_one.uns['cluster'], size=50, alpha=.1,min_edge_width=2, node_size_scale=1.5,save = path+'.paga.png',   add_outline=True, show=False)

#tempsample_one = sample_one.raw.to_sample_one()
#tempsample_one.var_names = sample_one.var_names
#sample_one.raw = tempsample_one
#scv.tl.score_genes_cell_cycle(sample_one)
#scv.pl.scatter(sample_one.raw, color_gradients=['S_score', 'G2M_score'], smooth=True, perc=[5, 95],save = file+'cycling.png',show=False)

#######################################################################################
