#################################################################
#
# date: April 12, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# DE-speudoBulk.R
# Usage : 
# 
# speudoBulk.R ~TODO
# 
# Description : 
#
# SpeudoBulk.R differetial Epresssion (DE) of scRNA-SEQ
#
#################################################################

library(plyr)
library(dplyr)
library(Seurat)
library(glue)
library(hdf5r)
library(scater) # Davis McCarthy Bioinformatics 2017 Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R
library(patchwork)
library(SingleCellExperiment)# Robert A. Amzequita 2019 Nature merthods Orchestrating single-cell analysis with Bioconductor
library(scran)
library(org.Hs.eg.db)
library(stringr)
library(biomaRt)


#######################################################################################################
###################################       Load h5 files        ########################################
#######################################################################################################
plotting = FALSE
set.seed(100)
# nCount_RNA = the number of UMIs per cell nUMI
# nFeature_RNA = the number of genes detected per cell nGene
# ------------------------------------------------------------------------------------------------------------
#cond1 = "OSI_TIPI_Vertes"
#cond2 = "OSI_TIPI_Rouges"

cond1 = "CTL_Vertes"
cond2 = "CTL_Rouges"

base.dir = "/data/villemin/data/toulouse/scRNAseqCells-2/CellRanger"

#cond1 = PC9_verte
#cond2 = PC9_rouge

#cond1 = "4006_verte"
#cond2 = "4006_rouge"

#base.dir = "/data/villemin/data/toulouse/scRNAseqCells-1/CellRanger"
# ------------------------------------------------------------------------------------------------------------
dir.create(glue("{base.dir}/DE"), showWarnings = F)

#sce.to.seurat <- readRDS(file = glue("{base.dir}/{cond1}_vs_{cond2}.rds"))
sce.to.seurat <- readRDS(file = glue("{base.dir}/CTL.test.rds"))


#######################################################################################################
###################################      DE speudo bulk           #####################################
#######################################################################################################


summed <- aggregateAcrossCells(seurat.to.sce, id = colData(merged)[,c("orig.ident")])
summed
                                                                                                          