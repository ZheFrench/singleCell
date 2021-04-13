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
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# Description : 
#
# SpeudoBulk.R differetial Epresssion (DE) of scRNA-SEQ
#
#################################################################
library(plyr)
library(dplyr)
library(tidyverse)

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
library(magrittr)
library(edgeR)
library(data.table)
library(Matrix.utils)

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


if (FALSE){
    
print(" Loading RDS.. ")
    
seurat.Object <- readRDS(file = glue("{base.dir}/{cond1}_vs_{cond2}.rds"))
print(head(seurat.Object[[]])) # data@meta.data

#######################################################################################################
###################################      DE speudo bulk           #####################################
#######################################################################################################
    
#Idents(cd14.mono) <- "stim"
    
# Switch seurat object to SingleCellExperiment (Nat methods 2020)
sce.Object <- as.SingleCellExperiment(seurat.Object)
sce.Object$ident <- sce.Object$orig.ident

length(sce.Object$ident)#5775
# Check the assays present
#print(head(colData(sce.Object)))

nrows.cond1 = dim(sce.Object[, sce.Object$ident == cond1])[2]
nrows.cond2 = dim(sce.Object[, sce.Object$ident == cond2])[2]
print(nrows.cond1)
print(nrows.cond2)

#print(length(sce.Object[, sce.Object$ident==cond2]))

random1 <- sample.int(3, nrows.cond1,replace = TRUE,prob=c(0.33,0.33,0.33))
random2 <- sample.int(3, nrows.cond2,replace = TRUE,prob=c(0.33,0.33,0.33))

head(random1)
head(random2)

sce.Object$subgroup <-"test"

# VERY LONG PROCESS
print ("set random1")
sce.Object[, sce.Object$ident==cond1]$subgroup <- random1
print("set random2")
sce.Object[, sce.Object$ident==cond2]$subgroup <- random2

# Save a seurat object... 
saveRDS(sce.Object,file = glue("{base.dir}/{cond1}_vs_{cond2}_subgroup.rds"))
    
}

sce.Object <- readRDS(file = glue("{base.dir}/{cond1}_vs_{cond2}_subgroup.rds"))

print("##### INITIAL ######")
print("coldData")
print(head(colData(sce.Object)))

print("rowData")
print(head(rowData(sce.Object)))

# remove lowly expressed gene
dim(sce.Object) #33694
sce.Object <- sce.Object[rowSums(counts(sce.Object) > 1) >= 150, ]
dim(sce.Object) #12247

sce.Object$sample_id    <- paste(sce.Object$orig.ident, "-", sce.Object$subgroup,sep = "")

# just in order to get the same names as tutorial
#sce.Object$group.id    <- sce.Object$orig.ident
# cluster.id  <- orig.ident equals group.id in tutorial

nk <- length( kids <- sce.Object$orig.ident[!duplicated(sce.Object$orig.ident)] )         
ns <- length( sids <- unique(sce.Object$sample_id) )

#print (nk)
#print (ns)

m <- match(sids, sce.Object$sample_id)
#print(m)
#3334 3333 3339    7    1    3
n_cells <- as.numeric(table(sce.Object$sample_id))
ei <- data.frame(colData(sce.Object)[m, ], n_cells, row.names = NULL) 
print("ei")
print(head(ei))

# Aggregate across cluster-sample groups
groups <- colData(sce.Object)[, c("sample_id")]
# split by cluster, transform & rename columns
matrix.rnaseq <- aggregate.Matrix(t(counts(sce.Object)), groupings = groups, fun = "sum") 

#http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/OSCABioc2019__OSCABioc2019/
#http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/muscWorkshop__vignette/
#https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html

#matrix.rnaseq <- split(matrix.rnaseq, unique(sce.Object$sample_id))  %>% lapply(function(u)  set_colnames(t(u), unname(sids))  )
                                                              
# construct SCE of pseudo-bulk counts
matrix.rnaseq.cse <- SingleCellExperiment(assays =list(counts=t(matrix.rnaseq))  )

matrix.rnaseq.cse
# Add QC metrics from scater
colData(matrix.rnaseq.cse) <- cbind(colData(matrix.rnaseq.cse),perCellQCMetrics(matrix.rnaseq.cse)) # sum / detected / total
rowData(matrix.rnaseq.cse) <- cbind(rowData(matrix.rnaseq.cse),perFeatureQCMetrics(matrix.rnaseq.cse)) # mean / detected

print("##### SPEUDO BULK ######")

print("coldData")
print(head(colData(matrix.rnaseq.cse)))

print("rowData")
print(head(rowData(matrix.rnaseq.cse)))

stop()




sce.filt <- logNormCounts(matrix.rnaseq.cse)

dec <- modelGeneVar(sce.filt, block = sce.filt$sample_id)
hvgs = getTopHVGs(dec, n = 2000)

sce.filt <- runPCA(sce.filt, use_coldata = TRUE, detect_outliers = TRUE)

png(file = glue("{base.dir}/plots/{cond1}_{cond2}_PCA.png"),width = 1500,height = 500)
plotReducedDim(sce.filt, use_dimred="PCA_coldata", colour_by = "ident")
dev.off()

png(file = glue("{base.dir}/plots/Dim_PCA.png"),width = 1500,height = 500)
vdl <- VizDimLoadings(object = sce.filt, dims = 1:3)
print(vdl)
dev.off()
# Creating up a DGEList object for use in edgeR:

#y <- DGEList(counts(sce.Object), samples=colData(sce.Object))
