import numpy as np
import pandas as pd
from anndata import AnnData
from scanpy import read
import scvelo as scv
import loompy

# First 
# velocyto2.R
# Then scVelo
# export HDF5_USE_FILE_LOCKING='FALSE'
# ssh -X

file1 = "/data/USERS/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_CTL_GRCh38/velocyto/Human_CTL_GRCh38.loom"
file2 = "/data/USERS/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_OSI/velocyto/Human_OSI.loom"
file3 = "/data/USERS/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_OSI_TIPI/velocyto/Human_OSI_TIPI.loom"
file4 = "/data/USERS/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_TIPI/velocyto/Human_TIPI.loom"

loompy.combine([file1,file2,file3,file4], key = "Accession",output_file = "PDX_ALL.loom")

exit(0)


file1 = "/data/USERS/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_rouge/velocyto/4006_rouge.loom"
file2 = "/data/USERS/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_verte/velocyto/4006_verte.loom"

#loompy.combine([file1,file2], key = "Accession",output_file = "4006_OSI.loom")

#file.combined = "/data/villemin/code/singleCell/python/4006_OSI.loom"

file3 = "/data/USERS/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Rouges/velocyto/CTL_Rouges.loom"
file4 = "/data/USERS/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Vertes/velocyto/CTL_Vertes.loom"

#loompy.combine([file3,file4], key = "Accession",output_file = "4006_CTL.loom")

#loompy.combine([file1,file2,file3,file4], key = "Accession",output_file = "4006_CTL_OSI.loom")


file5 = "/data/USERS/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Rouges/velocyto/OSI_TIPI_Rouges.loom"
file6 = "/data/USERS/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Vertes/velocyto/OSI_TIPI_Vertes.loom"

loompy.combine([file5,file6], key = "Accession",output_file = "4006_OSI_TIPI.loom")

loompy.combine([file1,file2,file3,file4,file5,file6], key = "Accession",output_file = "4006_ALL.loom")

