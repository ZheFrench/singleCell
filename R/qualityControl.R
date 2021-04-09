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
#################################################################
library(optparse)
option_list = list(
  make_option(c("-a", "--condition1"), type="character", default=NULL, help="Condition1 is the name of output directory of CellRanger for one condition.", metavar="CONDITION"),
  make_option(c("-b", "--condition2"), type="character", default=NULL, help="Condition1 is the name of output directory of CellRanger for another condition.", metavar="CONDITION"),
  make_option(c("-d", "--directory"), type="character", default=NULL, help="Directory to look for output directories of CellRanger. ", metavar="PATH2DIRECTORY"),
  make_option(c("-p", "--plot"), type="logical", default=TRUE, help="Plot things or not. ", metavar="TRUE|FALSE")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = 0);
opt <- arguments$options

print("> ARGS : ")
print(opt$directory)
print(opt$condition1)
print(opt$condition2)
print(opt$plot)

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
plotting = opt$plot
set.seed(100)
# nCount_RNA = the number of UMIs per cell nUMI
# nFeature_RNA = the number of genes detected per cell nGene
# ------------------------------------------------------------------------------------------------------------
#cond1 = "OSI_TIPI_Vertes"
#cond2 = "OSI_TIPI_Rouges"

cond1 = opt$condition1
cond2 = opt$condition2

base.dir = opt$directory

#cond1 = PC9_verte
#cond2 = PC9_rouge

#cond1 = "4006_verte"
#cond2 = "4006_rouge"

#base.dir = "/data/villemin/data/toulouse/scRNAseqCells-1/CellRanger"
# ------------------------------------------------------------------------------------------------------------
dir.create(glue("{base.dir}/plots"), showWarnings = F)

# Read 10x output
data_condition1 = Read10X_h5(glue("{base.dir}/{cond1}/outs/filtered_feature_bc_matrix.h5"))
data_condition2 = Read10X_h5(glue("{base.dir}/{cond2}/outs/filtered_feature_bc_matrix.h5"))

dim(x = data_condition1)
dim(x = data_condition2)

# Create individual Seurat objects
cond1.object = CreateSeuratObject(data_condition1, project = cond1)
cond2.object = CreateSeuratObject(data_condition2, project = cond2)

cond1.object$log10GenesPerUMI <- log10(cond1.object$nFeature_RNA) / log10(cond1.object$nCount_RNA)
cond2.object$log10GenesPerUMI <- log10(cond2.object$nFeature_RNA) / log10(cond2.object$nCount_RNA)


# Compute percent mito ratio
cond1.object$mitoRatio <- PercentageFeatureSet(object = cond1.object, pattern = "^MT-")
cond2.object$mitoRatio <- PercentageFeatureSet(object = cond2.object, pattern = "^MT-")

data <- merge(cond1.object, cond2.object, add.cell.ids=c(cond1,cond2))

#head(data@meta.data)
head(data@meta.data)

#Low library size: When cells are very degraded or absent from the library preparation, the number of reads sequenced from that library will be very low. It’s important to remove these cells from downstream analyses.

#Low number of expressed genes: A low number of expressed genes may be a result of poor-quality cells (e.g. dying, degraded, damaged, etc.), followed by high PCR amplification of the remaining RNA. Again, these cells should be removed from downstream analyses.

#High mitochondrial gene content: High concentrations of mitochondrial genes is often a result of damaged cells where the endogenous RNA escapes or degrades. As mitochondria has its own cell membranes, it is often the last DNA/RNA in damaged cells to degrade and hence occurs in high quantities during sequencing.

#Batch effect: Large scRNA-seq projects usually need to generate data across multiple batches due to logistical constraints. However, the processing of different batches is often subject to variation, e.g., changes in operator, differences in reagent quality and concentration, the sequencing machine used, etc. This results in systematic differences in the observed expression in cells from different batches, which we refer to as “batch effects”. Batch effects are problematic as they can be major drivers of variation in the data, masking the relevant biological differences and complicating interpretation of the results.
# Filter out low quality reads using selected thresholds - these will change with experiment

# -------------------------
# Create metadata dataframe
metadata <- data@meta.data


#######################################################################################################
###################################       Seurat QC        ############################################
#######################################################################################################
    

if (plotting){

# -------------------------
# Visualize the number UMIs/transcripts per cell
density.plot.nCount_RNA <- metadata %>% 
  ggplot(aes(color=orig.ident, x = nCount_RNA, fill = orig.ident)) + 
   geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
  # + geom_vline(xintercept = 500)

density.plot.nFeature_RNA <- metadata %>% 
  ggplot(aes(color=orig.ident, x = nFeature_RNA, fill = orig.ident)) + 
   geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density")

density.plot.log10GenesPerUMI <- metadata %>% 
  ggplot(aes(color=orig.ident, x = log10GenesPerUMI, fill = orig.ident)) + 
   geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 

density.plot.mitoRatio <- metadata %>% 
  ggplot(aes(color=orig.ident, x = mitoRatio, fill = orig.ident)) + 
   geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
 #+ geom_vline(xintercept = 500)

png(file = glue("{base.dir}/plots/{cond1}_vs_{cond2}_RNA_density.png"),width = 1500,height = 500)
print(density.plot.nCount_RNA + density.plot.nFeature_RNA  + density.plot.log10GenesPerUMI  + density.plot.mitoRatio )
dev.off()
    
# -------------------------

print("::: QC PLOT - Seurat :::")

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_violonQC.png"),width = 1500,height = 500)
print(VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio","log10GenesPerUMI"), ncol = 4)  + theme(  axis.text = element_text( size = 14)))
dev.off()

# Visualize the distribution of cell cycle markers across
png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_cellcycle_ridgleplot.png"),width = 1500,height = 500)
print(RidgePlot(data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2))
dev.off()


}



data <- NormalizeData(data)

# Add Cell Cycle (need to be normalized before)
#https://github.com/satijalab/seurat/issues/1679
data <- CellCycleScoring( object = data,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)#,  set.ident = TRUE
metadata <- data@meta.data
print(head(metadata))

if (plotting){
    
png(file = glue("{base.dir}/plots/{cond1}_vs_{cond2}_cellcycle_score.png"),width = 1500,height = 500)
print(VlnPlot(data, features = c("S.Score","G2M.Score")))
dev.off()
#https://satijalab.org/seurat/articles/cell_cycle_vignette.html#assign-cell-cycle-scores
   
histo.plot.phase <- metadata %>% 
ggplot(aes(x=orig.ident, fill=Phase)) + 
geom_bar(alpha = 0.2,colour="black") +  
ggtitle("NCells")
    
png(file = glue("{base.dir}/plots/{cond1}_vs_{cond2}_phase_histo.png"),width = 1500,height = 500)
print(histo.plot.phase + theme (axis.text = element_text(size = 14 )))
dev.off()  
    
}

#######################################################################################################
###################################       Normalize            ########################################
#######################################################################################################
print ("::: Normalize :::")

data <- FindVariableFeatures(data, selection.method = "vst") #nFeatures=2000
data <- ScaleData(data, features = rownames(data))
     
#######################################################################################################
###################################       PCA            ########################################
#######################################################################################################

# Perform PCA

if (FALSE) { 
data <- JackStraw(data, num.replicate = 100)
data <- ScoreJackStraw(data, dims = 1:50)

png(file = glue("{base.dir}/plots/{cond1}_vs_{cond2}_dimension_pca.png"),width = 1500,height = 500)
j <- JackStrawPlot(data, dims = 1:50)
e <- ElbowPlot(data.norm)
print(e+j)
dev.off() 
}

if (plotting){

data <- RunPCA(data, features = VariableFeatures(data))

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_genesBigVariance_pca.png"),width = 1500,height = 500)
print(DimPlot(data,shape.by="orig.ident",group.by = "Phase",raster = FALSE,pt.size = 2))
dev.off() 
    
data <- RunPCA(data, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
    
png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_genesCellCycle_pca.png"),width = 1500,height = 500)
print(DimPlot(data,shape.by="orig.ident",group.by = "Phase",raster = FALSE,pt.size = 2) )   
dev.off() 

}

print("Mydata after Normalization & PCA :")# data@meta.data
print(head(data[[]]))# data@meta.data


#######################################################################################################
###################################       SingleCellExperiment QC   (scater)     ######################
#######################################################################################################


# Switch seurat object to SingleCellExperiment (Nat methods 2020)
seurat.to.sce <- as.SingleCellExperiment(data)

# Check the assays present
#assays(seurat.to.sce)

#head(seurat.to.sce)
    
# Add QC metrics from scater
colData(seurat.to.sce) <- cbind(colData(seurat.to.sce),perCellQCMetrics(seurat.to.sce)) # sum / detected / total
# total is the total UMIs per cell. nUMi, nCount
# detected is the number of detected genes per cell . nFeature
# sum ?

rowData(seurat.to.sce) <- cbind(rowData(seurat.to.sce),perFeatureQCMetrics(seurat.to.sce)) # mean / detected
# mean is the mean UMI/count per cell.
# detected is the number of umi/count for this feature.

print(head(colData(seurat.to.sce)))
print(head(rowData(seurat.to.sce)))

if (plotting){
    
print ("::: QC PLOT - scater :::")

# Plot an overview of expression for each cell
png(file = glue("{base.dir}/plots/{cond1}_vs_{cond2}_n_cells_by_counts.png"),width = 1000,height = 500)
print(plotRowData(seurat.to.sce, y="detected", x="mean") +  scale_x_log10())
dev.off()

# Plot the percentage of counts accounted for by the top most highly expressed features across the dataset.
# Each row on the plot corresponds to a feature and is sorted by average expression
# These ticks can be coloured according to cell-level metadata, as specified by colour_cells_by

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_scater_mostexpresed.png"),width = 1600,height = 750)

p1 <- plotHighestExprs(seurat.to.sce[, seurat.to.sce$ident==cond1], exprs_values = "counts",colour_cells_by="detected")   +  theme(  axis.text = element_text( size = 14), legend.text = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_blank(),legend.title = element_text( size = 14 )  ) + 
 xlab("Percentage of counts") + ylab("Genes") +labs(fill = "Counts per Cell")

p2 <- plotHighestExprs(seurat.to.sce[, seurat.to.sce$ident==cond2], exprs_values = "counts",colour_cells_by="detected")  + theme( axis.text = element_text( size = 14), legend.text = element_text(size=14), axis.title.x = element_text(size=14),axis.title.y = element_blank(),legend.title = element_text( size = 14 ) ) +
xlab("Percentage of counts") + ylab("Genes") + guides(fill=guide_legend(title="Counts per Cell"))
      
print(p1 + p2 + plot_layout(ncol=2) )
      
dev.off()

}

#theme(base_size = 22)



#######################################################################################################
###################################      Filtering            #########################################
#######################################################################################################
print ("::: Filtering (TODO) :::")

#filtered_seurat.cond1 <- subset(x = cond1.object,subset= (nUMI >= 9000) &(nGene >= 2800) &(log10GenesPerUMI >= 0.8) & (mitoRatio < 0.20))
# 4006_red
#filtered_seurat.cond1 <- subset(x = cond1.object,subset= (nUMI >= 9000) & (nGene >= 2800) & (log10GenesPerUMI >= 0.8) & (mitoRatio < 0.20))
# 4006_green
#filtered_seurat.cond2 <- subset(x = cond2.object, subset= (nUMI >= 9000) & (nGene >= 2800) & (log10GenesPerUMI >= 0.8) & (mitoRatio < 0.20))


# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
#counts <- GetAssayData(object = filtered_seurat, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
#nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
#keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
#filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
#filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

# Merge two seurat objects by its idents
#TwoConds.object <- merge(filtered_seurat.cond1, filtered_seurat.cond2)
#filtered_seurat.cond1 = TwoConds.object


# To remove specific genes from the seurat object ==================
#genes_to_remove <- fread("H4006_removed_genes.txt",data.table=F) # Given list of genes to be removed.

#genes_to_remove = as.vector(unlist(genes_to_remove))
#counts <- GetAssayData(TwoConds.object, assay = "RNA")
#counts <- counts[-(which(rownames(counts) %in% genes_to_remove)),]
#Patient.object <- subset(TwoConds.object, features = rownames(counts))
#filtered_seurat.cond1 = Patient.object

# Filter out low quality reads using selected thresholds - these will change with experiment
#keep <- metrics %>%
#  dplyr::filter(nUMI > 500 , 
#               nGene > 500,
#                log10GenesPerUMI > 0.8,
#               mitoRatio < 0.1,
#                ) %>% 
 # pull(cells)

# Subset the cells to only include those that meet the thresholds specified
#se_c <- se[ ,keep]

# Save subset to new metrics variable
#metrics_clean <- colData(se_c) %>%
# as.data.frame()

# Save cleaned single-cell experimnet as .RData to load at any time
#saveRDS(se_c, file = "data/se_filtered.rds")

#######################################################################################################
###################################       Calculate cell-cycle scores            ######################
#######################################################################################################
print ("::: cell-cycle scores with cyclone from   :::")

if (TRUE){
    
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))


organism = "hsapiens_gene_ensembl"
host = "mar2016.archive.ensembl.org" # ensembl 84 ... There is no 83 ...
symbol_description = 'hgnc_symbol'

# BIOMART OBJECT
edb = useMart("ENSEMBL_MART_ENSEMBL", dataset=organism,host=host)

# Retrieve gene infos
gene_infos = getBM(attributes=c('ensembl_gene_id',symbol_description,'gene_biotype'),values=rownames(seurat.to.sce),filters=symbol_description,mart=edb)

# Output had duplicate lines
gene_infos_without_dup <- gene_infos[!duplicated(gene_infos$hgnc_symbol),]

df.genes <- data.frame(hgnc_symbol=rownames(seurat.to.sce))

res_annotated <- join(df.genes , gene_infos_without_dup, by = symbol_description, type='left', match='all') #first

rowData(seurat.to.sce)$feature_ensembl <-res_annotated$ensembl_gene_id

ensembl <- res_annotated$ensembl_gene_id
# Use only genes related to biological process cell cycle to speed up
# https://www.ebi.ac.uk/QuickGO/term/GO:0007049 = cell cycle (BP,Biological Process)
GOs <- na.omit(select(org.Hs.eg.db, keys = na.omit(ensembl), keytype = "ENSEMBL", column = "GO"))
GOs <- GOs[GOs$GO == "GO:0007049", "ENSEMBL"]
hs.pairs <- lapply(hs.pairs, function(x) { x[rowSums(apply(x, 2, function(i) i %in% GOs)) >= 1, ]}) # don't know why
str(hs.pairs)
#cc.ensembl <- ensembl[ensembl %in% GOs] # 250 fast but less accurate
cc.ensembl <- ensembl[ ensembl %in% unique(unlist(hs.pairs))] # 1189

print("hs.pairs")
print(length(cc.ensembl ))
# define feature names in feature_symbol column
rowData(seurat.to.sce)$feature_symbol <- rownames(seurat.to.sce)

# remove features with duplicated names
#seurat.to.sce <- seurat.to.sce[!duplicated(rowData(seurat.to.sce)$feature_symbol), ]

# Using Ensembl IDs to match up with the annotation in 'hs.pairs'.
assignments <- cyclone(seurat.to.sce[rowData(seurat.to.sce)$feature_ensembl %in% cc.ensembl, ], hs.pairs, gene.names =  ensembl[ensembl %in% cc.ensembl])
                                                    
seurat.to.sce$phases.cyclone    <- assignments$phases
seurat.to.sce$G1.score.cyclone  <- assignments$scores$G1
seurat.to.sce$G2M.score.cyclone <- assignments$scores$G2M
seurat.to.sce$S.score.cyclone   <- assignments$scores$S
              
#table(assignments$phases, (seurat.to.sce))

png(file = glue("{base.dir}/plots/{cond1}_vs_{cond2}_cellcyle_score_cyclone.png"),width = 1000,height = 500)
                                                    
p2 <-plotColData(seurat.to.sce, y = "G2M.score.cyclone", x = "orig.ident", colour_by = "orig.ident")
p3 <-plotColData(seurat.to.sce,y = "G1.score.cyclone", x = "orig.ident", colour_by = "orig.ident")
p4 <- plotColData(seurat.to.sce, y = "S.score.cyclone", x = "orig.ident", colour_by = "orig.ident")
                                                     
print( p2 + p3 + p4 )
 #https://satijalab.org/seurat/articles/cell_cycle_vignette.html#assign-cell-cycle-scores
           
dev.off()      
                                                 
histo.plot.phase <- as.data.frame(colData(seurat.to.sce)) %>% 
ggplot(aes(x=orig.ident, fill=phases.cyclone)) + 
geom_bar(alpha = 0.2,colour="black") +  
ggtitle("NCells")
    
png(file = glue("{base.dir}/plots/{cond1}_vs_{cond2}_phase_histo_cyclone.png"),width = 1500,height = 500)
print(histo.plot.phase + theme (axis.text = element_text(size = 14 )))
dev.off()  
}
                                                     
# ----------------------------------------------------------------------------------
# Go back to seurat - It doesn't bring  colDat and rowData
#sce.to.seurat <- as.Seurat(seurat.to.sce)# Doesn"t work 
                                                     
# https://www.biostars.org/p/9464198/#9464336
# # I keep the meta.data from seurat , didn't succed to get back rowData so I reprocess seurat object.
sce.to.seurat <- CreateSeuratObject(counts = counts(seurat.to.sce), meta.data = as.data.frame(colData(seurat.to.sce)))
#print(head(data[[]]))# data@meta.data
 
# Redo that because we loose stuffs passing from sce to seurat
sce.to.seurat <- NormalizeData(sce.to.seurat)
sce.to.seurat <- FindVariableFeatures(sce.to.seurat, selection.method = "vst") #nFeatures=2000
sce.to.seurat <- ScaleData(sce.to.seurat, features = rownames(sce.to.seurat))
# ----------------------------------------------------------------------------------
   
sce.to.seurat <- RunPCA(sce.to.seurat, features = VariableFeatures(sce.to.seurat))

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_genesBigVariance_pca_cyclone.png"),width = 1500,height = 500)
print(DimPlot(sce.to.seurat,shape.by="orig.ident",group.by = "phases.cyclone",raster = FALSE,pt.size = 2))
dev.off() 
    
sce.to.seurat <- RunPCA(sce.to.seurat, features = c(cc.genes$s.genes, cc.genes$g2m.genes))
    
png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_genesCellCycle_pca_cylone.png"),width = 1500,height = 500)
print(DimPlot(sce.to.seurat,shape.by="orig.ident",group.by = "phases.cyclone",raster = FALSE,pt.size = 2) )   
dev.off()                                                     
                                                  
# Save a seurat object... 
saveRDS(sce.to.seurat,file = glue("{base.dir}/{cond1}_vs_{cond2}.rds"))