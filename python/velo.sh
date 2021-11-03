#!/bin/bash
export HDF5_USE_FILE_LOCKING='FALSE'

##########################################
######			Velocity			######
##########################################

################# PDX ################
python /data/villemin/code/singleCell/python/loom_filter.py -l  /data/villemin/code/singleCell/python/PDX_ALL.loom -t /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/clusters/Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI_cellID_obs.csv

python /data/villemin/code/singleCell/python/velocyto_1.py -f  /data/villemin/code/singleCell/python/PDX_ALL.loom   -c Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI -d "/data/villemin/data/toulouse/scRNAseqPDX/CellRanger/clusters/"

#Rscript /data/villemin/code/singleCell/R/plotNestedCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI
Rscript /data/villemin/code/singleCell/R/ssGsea.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI

exit 0
################# 4006_OSI_TIPI ################

python /data/villemin/code/singleCell/python/loom_filter.py -l  /data/villemin/code/singleCell/python/4006_ALL.loom -t /data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/4006_rouge.4006_verte.CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes_cellID_obs.csv

python /data/villemin/code/singleCell/python/velocyto_1.py -f  /data/villemin/code/singleCell/python/4006_ALL.loom   -c 4006_rouge.4006_verte.CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes -d "/data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/"

exit 0

################# 4006_CTL ################

python /data/villemin/code/singleCell/python/loom_filter.py -l  /data/villemin/code/singleCell/python/4006_CTL.loom -t /data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/CTL_Rouges.CTL_Vertes_cellID_obs.csv

python /data/villemin/code/singleCell/python/velocyto_1.py -f  /data/villemin/code/singleCell/python/4006_CTL.loom   -c CTL_Rouges.CTL_Vertes -d "/data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/"

exit 0
################# 4006_OSI ################

python /data/villemin/code/singleCell/python/loom_filter.py -l  /data/villemin/code/singleCell/python/4006_OSI.loom -t /data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/4006_rouge.4006_verte_cellID_obs.csv

python /data/villemin/code/singleCell/python/velocyto_1.py -f  /data/villemin/code/singleCell/python/4006_OSI.loom   -c 4006_rouge.4006_verte -d "/data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/"

python /data/villemin/code/singleCell/python/velocyto_2.py -f  /data/villemin/code/singleCell/python/4006_OSI_filtered.h5ad

# Super long 
Rscript /data/villemin/code/singleCell/R/velo.R -f /data/villemin/code/singleCell/python/4006_OSI_filtered.loom

exit 0
################# 4006_CTL_OSI ################

python /data/villemin/code/singleCell/python/loom_filter.py -l  /data/villemin/code/singleCell/python/4006_CTL_OSI.loom 
-t /data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/4006_rouge.4006_verte.CTL_Rouges.CTL_Vertes_cellID_obs.csv


python /data/villemin/code/singleCell/python/velocyto_1.py -f  /data/villemin/code/singleCell/python/4006_CTL_OSI.loom   -c 4006_rouge.4006_verte.CTL_Rouges.CTL_Vertes -d "/data/villemin/data/toulouse/scRNA-ALL/CellRanger/clusters/"

python /data/villemin/code/singleCell/python/velocyto_2.py -f  /data/villemin/code/singleCell/python/4006_CTL_OSI_filtered.h5ad

# Super long 
Rscript /data/villemin/code/singleCell/R/velo.R -f /data/villemin/code/singleCell/python/4006_CTL_OSI_filtered.loom

##########################################
######			Cluster			######
##########################################
exit
# PREPROCESS CLUSTER / UMAP 
Rscript /data/villemin/code/singleCell/R/plotNestedCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger CTL_Rouges CTL_Vertes
Rscript /data/villemin/code/singleCell/R/plotNestedCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger 4006_rouge 4006_verte
Rscript /data/villemin/code/singleCell/R/plotNestedCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger 4006_rouge 4006_verte CTL_Rouges CTL_Vertes
 Rscript /data/villemin/code/singleCell/R/plotNestedCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger 4006_rouge 4006_verte CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes

# DIFF EXP
Rscript /data/villemin/code/singleCell/R/diffCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger 4006_rouge 4006_verte
#Rscript /data/villemin/code/singleCell/R/diffCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger 4006_rouge 4006_verte CTL_Rouges CTL_Vertes
#Rscript /data/villemin/code/singleCell/R/diffCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger CTL_Rouges CTL_Vertes

# PDX 
Rscript /data/villemin/code/singleCell/R/plotNestedCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI
#Rscript /data/villemin/code/singleCell/R/diffCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI
