import numpy as np
import pandas as pd
from anndata import AnnData
from scanpy import read
import scvelo as scv
import loompy
from pathlib import Path
import argparse
import os


parser = argparse.ArgumentParser()
parser.add_argument("-l", "--loom", action="store", help="file loom", required=True, type=str, dest='loom'),
parser.add_argument("-t", "--threshold", action="store", help="file ids cell", required=True, type=str, dest='threshold')

parameters = parser.parse_args()

print(parameters.loom)
print(parameters.threshold)

filename = Path(parameters.loom)
filename_wo_ext = str(filename.with_suffix(''))

file = os.path.basename(parameters.loom)
print(file)
print(filename_wo_ext)

# First 
# velocyto2.R
# Then scVelo
# export HDF5_USE_FILE_LOCKING='FALSE'
# ssh -X
print("sample_obs")
sample_obs    = pd.read_csv(parameters.threshold)
sample_obs.rename(columns={'Cells(seurat.Object)':'CellID'}, inplace=True )

sample_obs['CellID'] = sample_obs['CellID'].str[::-1]
sample_obs['CellID'] = sample_obs['CellID'].str.replace('_',':',1)
sample_obs['CellID'] = sample_obs['CellID'].str[::-1]
sample_obs['CellID'] = sample_obs['CellID'].str.replace('-1','x')

print(sample_obs.head())

outputname = filename_wo_ext+"_filtered.loom"

with loompy.connect(parameters.loom) as ds:
    view = ds.view[:, np.isin(ds.ca.CellID, sample_obs['CellID'])]
    
loompy.create(outputname, view.layers, view.ra, view.ca)