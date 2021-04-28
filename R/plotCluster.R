#################################################################
#
# date: April 02, 2021
# platform: Ubuntu 10.04
# R.version : 4.0.3
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# plotClusters.R
# Usage : 
# 
# plotClusters directory condition1 condition2
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

conditions.all <- strsplit(input,".",fixed = TRUE)

subset.string      <- ""


if (length(as.list(args)) == 0){
    subset.string <-   opt$input 
    conditions.to.test <- strsplit(input,".",fixed = TRUE)
}  else { subset.string <- paste0(conditions.to.test, collapse = ".")  }

  
combined <- readRDS(file = glue("{base.dir}/{subset.string}_integrated.rds"))

head(combined@meta.data)
#str(combined)

# Visualization[, combined$orig.ident %in% conditions.to.test]
p1 <- DimPlot(combined, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "umap", label = TRUE, repel = TRUE)
p3 <- DimPlot(combined, reduction = "umap", group.by = "Phase")

png(file = glue("{base.dir}/clusters/{subset.string}_umap_cluster.png"),width = 1800,height = 600)
p1 | p2 | p3
dev.off()

combined <- RunTSNE(object = combined)

# Visualization[, combined$orig.ident %in% conditions.to.test]
p1 <- DimPlot(combined, reduction = "tsne", group.by = "orig.ident")
p2 <- DimPlot(combined, reduction = "tsne", label = TRUE, repel = TRUE)
p3 <- DimPlot(combined, reduction = "tsne", group.by = "Phase")

png(file = glue("{base.dir}/clusters/{subset.string}_tnse_cluster.png"),width = 1800,height = 600)
p1 | p2 | p3
dev.off()


# Cells 
#rouges <- c ("CTL_Rouges","4006_rouge","OSI_TIPI_Rouges")
#vertes <- c ("CTL_Vertes","4006_verte","OSI_TIPI_Vertes")

#head(combined@meta.data)

# For performing differential expression after integration, we switch back to the original data
DefaultAssay(combined) <- "RNA"

print(length(conditions.to.test))
if (length(conditions.to.test) == 2 ) {

# order of comp
subset.string.corrected <- gsub(".","_", subset.string,fixed = TRUE)
print(subset.string.corrected)
print(conditions.to.test)

subset.string.corrected <- gsub(".","_", subset.string,fixed = TRUE)
 
if (  file.exists(glue("{base.dir}/DE/{subset.string.corrected}-differential-up.tsv")) != TRUE  ) {
    print(rev(conditions.to.test))
    subset.string.corrected <- gsub(".","_", paste0(rev(conditions.to.test), collapse = "."),fixed = TRUE)

}
    
genes.diff.up   <- fread(glue("{base.dir}/DE/{subset.string.corrected}-differential-up.tsv"),data.table = F)
symbol.genes.up <- head( genes.diff.up[order(genes.diff.up$FDR),]$genes, 12)

genes.diff.down <- fread(glue("{base.dir}/DE/{subset.string.corrected}-differential-down.tsv"),data.table = F)
symbol.genes.down <-head( genes.diff.down[order(genes.diff.down$FDR),]$genes, 12)
    
print("Up")
print(symbol.genes.up)
print("Down")
print(symbol.genes.down)
typeof(symbol.genes.down)
    #combined[, combined$orig.ident %in% conditions.to.test]
png(file = glue("{base.dir}/clusters/{subset.string}_highlight_GenesDown.png"),width = 1200,height = 600)
print(FeaturePlot(combined, features = symbol.genes.down) )#, min.cutoff = "q9"
dev.off()
    
png(file = glue("{base.dir}/clusters/{subset.string}_highlight_GenesUp.png"),width = 1200,height = 600)
print(FeaturePlot(combined, features = symbol.genes.up) )#, min.cutoff = "q9"
dev.off()    
    
subsetGenes <- c("LMNB1", "CENPF", "CENPE", "BIRC5", "FOSL1", "ANLN", "FABP5", "MT2A", "AURKB", "MYOZ1", "CTGF", "ANKRD1", "BMF", "MKI67")
    
png(file = glue("{base.dir}/clusters/{subset.string}_subset.png"),width = 1200,height = 600)
print(FeaturePlot(combined, features = subsetGenes) )#, min.cutoff = "q9"
dev.off()
    

    
}
#markers.to.plot <- c("BIRC5", "FOSL1", "ANLN", "FABP5", "MT2A", "AURKB", "MYOZ1", "CTGF", "ANKRD1", "BMF")
#png(file = glue("{base.dir}/clusters/{input}_v4.png"),width = 1500,height = 750)
#DotPlot(combined, features = markers.to.plot, cols = c("blue", "red","orange", "green"), dot.scale = 8, split.by = "orig.ident") + RotatedAxis()
#dev.off()

# Heatmaps
#DoHeatmap(object = pbmc, features = heatmap_markers)

#Identify differential expressed genes across conditions (will use the edgeR approach)

#https://satijalab.org/seurat/reference/averageexpression

# Visualise....specific gene by Figure
#FeaturePlot(immune.combined, features = c("CD3D", "GNLY", "IFI6"), split.by = "stim", max.cutoff = 3, cols = c("grey", "red"))

#FindMarkers function will be used to find gene diff regulated between cells and clusters...

# Barplot...
#plots <- VlnPlot(immune.combined, features = c("LYZ", "ISG15", "CXCL10"), split.by = "stim", group.by = "celltype",pt.size = 0, combine = FALSE)
#wrap_plots(plots = plots, ncol = 1)
