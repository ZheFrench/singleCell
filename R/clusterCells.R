#################################################################
#
# date: April 02, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# qualityControl.R
# Usage : 
# 
# qualityControl.R directory condition1 condition2
# 
# Description : 
#
# Some quality controls for scRNAseq using Seurat, Scater, Scran etc etc....
#
# Hard Coded Options : 
#
# min.cells (default= 3 ) Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
# min.features (default= 200 ) Include cells where at least this many features are detected.
#################################################################
library(optparse)
option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, help="Directory to look for output directories of CellRanger. ", metavar="PATH2DIRECTORY"),
  make_option(c("-i", "--input"), type="character",  help="Input to rds ", metavar="PATH2FILE"),
  make_option(c("-s", "--subset"), default ="All" , type="character",  help="Condition to integrate together.", metavar="LISTCONDITIONS")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args


print("> OPTS : ")
print(opt$directory)
print(opt$input)
print("> ARGS : ")
print(args)


suppressPackageStartupMessages(library(plyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(glue))
suppressPackageStartupMessages(library(hdf5r))
suppressPackageStartupMessages(library(scater)) # Davis McCarthy Bioinformatics 2017 Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data 
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(SingleCellExperiment))# Robert A. Amzequita 2019 Nature merthods Orchestrating single-cell analysis with Bioconductor
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(metap))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tricycle))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))

#http://barc.wi.mit.edu/education/hot_topics/scRNAseq_2020/SingleCell_Seurat_2020.html

#######################################################################################################
#######################################################################################################
###################################      FUNCTIONS             ########################################
#######################################################################################################
#######################################################################################################
'%ni%' <- Negate('%in%')

#######################################################################################################
#######################################################################################################
###################################      MAIN             #############################################
#######################################################################################################
#######################################################################################################   

set.seed(100)

base.dir = opt$directory
input    = opt$input 
                                            
dir.create(glue("{base.dir}/clusters"), showWarnings = F)


conditions.to.test    <- as.list(args)
print(conditions.to.test)

subset.string      <- ""


if (length(as.list(args)) == 0){
    subset.string <-   opt$input 
    conditions.to.test <- strsplit(input,".",fixed = TRUE)
}  else { subset.string <- paste0(conditions.to.test, collapse=".")  }


seurat.Object <- readRDS(file = glue("{base.dir}/{input}.rds"))

list.objects.by.condition <- SplitObject(seurat.Object, split.by = "orig.ident")


index <-NULL
for (i in 1:length(list.objects.by.condition)){
     if ( names(list.objects.by.condition)[i] %ni% conditions.to.test){index <- rbind(index,i) }  
}
 
list.objects.by.condition[index] <- NULL      # Remove multiple list elements

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list.objects.by.condition  )

anchors <- FindIntegrationAnchors(object.list = list.objects.by.condition, anchor.features = features)

# this command creates an 'integrated' data assay
combined <- IntegrateData(anchorset = anchors)
                    
# specify that we will perform downstream analysis on the corrected data note that the original
# Unmodified data still resides in the 'RNA' assay
DefaultAssay(combined) <- "integrated"       
# Run the standard workflow for visualization and clustering
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunTSNE(object = combined)

saveRDS(combined,file = glue("{base.dir}/{subset.string}_integrated.rds"))
    
