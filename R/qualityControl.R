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
#
#################################################################

library(optparse)
option_list = list(
  make_option(c("-d", "--directory"), type="character", default=NULL, help="Directory to look for output directories of CellRanger. ", metavar="PATH2DIRECTORY"),
  make_option(c("-f", "--plot"), type="logical", default=TRUE, help="Plot things or not. ", metavar="TRUE|FALSE")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args


print("> OPTS : ")
print(opt$directory)
print(opt$plot)
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

#http://barc.wi.mit.edu/education/hot_topics/scRNAseq_2020/SingleCell_Seurat_2020.html

#######################################################################################################
#######################################################################################################
###################################      FUNCTIONS             ########################################
#######################################################################################################
#######################################################################################################

#######################################################################################################
###################################       Calculate cell-cycle scores            ######################
#######################################################################################################

cell.cycle <- function  (data,status,condition,seurat.to.sce){

print ("::: cell-cycle scores with cyclone from   :::")

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
p2 <-plotColData(seurat.to.sce, y = "G2M.score.cyclone", x = "orig.ident", colour_by = "orig.ident")
p3 <-plotColData(seurat.to.sce,y = "G1.score.cyclone", x = "orig.ident", colour_by = "orig.ident")
p4 <- plotColData(seurat.to.sce, y = "S.score.cyclone", x = "orig.ident", colour_by = "orig.ident")
     
png(file = glue("{base.dir}/plots/{status}_{condition}_cellcyle_score_cyclone.png"),width = 1000,height = 500)
print( p2 + p3 + p4 ) #https://satijalab.org/seurat/articles/cell_cycle_vignette.html#assign-cell-cycle-scores
dev.off()      
                                                 
histo.plot.phase <- as.data.frame(colData(seurat.to.sce)) %>% 
ggplot(aes(x = orig.ident, fill = phases.cyclone)) + 
geom_bar(alpha = 0.2,colour="black") +  
ggtitle("NCells")
    
png(file = glue("{base.dir}/plots/{status}_{condition}_phase_histo_cyclone.png"),width = 1500,height = 500)
print(histo.plot.phase + theme (axis.text = element_text(size = 14 )))
dev.off()  

                                                     
print("Normalize before calling Seurat CellCycleScoring")
data <- NormalizeData(data)
dim(x = data)

#cc.genes.updated.2019 - Version more recent but we are based on hold annotation so....

# Add Cell Cycle (need to be normalized before)
#https://github.com/satijalab/seurat/issues/1679
#https://www.bioinformatics.babraham.ac.uk/training/10XRNASeq/seurat_workflow.html - Look cool plot
data <- CellCycleScoring( object = data,g2m.features = cc.genes$g2m.genes,s.features = cc.genes$s.genes)#,  set.ident = TRUE

png(file = glue("{base.dir}/plots/{status}_{condition}_cellcycle_score.png"),width = 1500,height = 500)
print(VlnPlot(data, features = c("S.Score","G2M.Score")))
dev.off()
#https://satijalab.org/seurat/articles/cell_cycle_vignette.html#assign-cell-cycle-scores
   
histo.plot.phase <- data@meta.data %>% 
ggplot(aes(x=orig.ident, fill=Phase)) + 
geom_bar(alpha = 0.2,colour="black") +  
ggtitle("NCells")
    
png(file = glue("{base.dir}/plots/{status}_{condition}_phase_histo.png"),width = 1500,height = 500)
print(histo.plot.phase + theme (axis.text = element_text(size = 14 )))
dev.off()  

return (data)
}
                                                     
#############################################################################################################
###################################   Control.quality    ####################################################
#############################################################################################################
                                                                 
#Low library size: When cells are very degraded or absent from the library preparation, the number of reads sequenced from that library will be very low. It’s important to remove these cells from downstream analyses.

#Low number of expressed genes: A low number of expressed genes may be a result of poor-quality cells (e.g. dying, degraded, damaged, etc.), followed by high PCR amplification of the remaining RNA. Again, these cells should be removed from downstream analyses.

#High mitochondrial gene content: High concentrations of mitochondrial genes is often a result of damaged cells where the endogenous RNA escapes or degrades. As mitochondria has its own cell membranes, it is often the last DNA/RNA in damaged cells to degrade and hence occurs in high quantities during sequencing.

#Batch effect: Large scRNA-seq projects usually need to generate data across multiple batches due to logistical constraints. However, the processing of different batches is often subject to variation, e.g., changes in operator, differences in reagent quality and concentration, the sequencing machine used, etc. This results in systematic differences in the observed expression in cells from different batches, which we refer to as “batch effects”. Batch effects are problematic as they can be major drivers of variation in the data, masking the relevant biological differences and complicating interpretation of the results.
# Filter out low quality reads using selected thresholds - these will change with experiment


control.quality <- function  (data,status,condition){

    
###################################      GGPLOT QC        ###################################################

metadata <- data@meta.data
    
# Visualize Cells ordered by % ribosomal counts | % mitochondrial counts
ordered.gene.ribo <- metadata[order(metadata$riboRatio),] %>% ggplot(aes( x = seq_along(riboRatio), y=riboRatio, color = orig.ident)) + geom_point()  + facet_grid(. ~ orig.ident) +
  xlab("") +
  ylab("Ncounts") +
  ggtitle("Cells ordered by % ribosomal counts")

ordered.gene.mito <- metadata[order(metadata$mitoRatio),]  %>% ggplot(aes( x = seq_along(mitoRatio), y=mitoRatio, color = orig.ident)) + geom_point()  + facet_grid(. ~ orig.ident) +
  xlab("") +
  ylab("Ncounts") +
  ggtitle("Cells ordered by % mitochondrial counts")

png(file = glue("{base.dir}/plots/{status}_{condition}_mito_ribo_ordered_by_counts.png"),width = 1500,height = 500)
print(ordered.gene.ribo | ordered.gene.mito)
dev.off()

counts.cell <- as.data.frame(table(metadata$orig.ident))
colnames(counts.cell) <- c('Condition','Freq')


# Visualize the number of cell counts per cell
ccounts <- counts.cell %>% 
 ggplot(aes(x = Condition,y = Freq, fill = factor(Condition) )  ) + 
 geom_bar(alpha = 0.4 ,position="dodge",stat="identity",color="gray") +   coord_flip() +
 ggtitle("NCells")

png(file = glue("{base.dir}/plots/{status}_{condition}_totalCells_per_condition.png"),width = 1500,height = 500)
print(ccounts)
dev.off()


# Counts per cell
counts.per.cell <-  ggplot(metadata, aes(x = log10(as.numeric(nCount_RNA)+1)) )  +
       geom_histogram(aes(color = orig.ident, fill = orig.ident),  position = "identity", bins = 30, alpha = 0.4 )  + facet_grid(. ~ orig.ident) +
         xlab("") +
  ylab("log10(counts+1)") +
  ggtitle("Histogram :: Counts(UMI) per cell")

# Genes per cell
genes.per.cell <-  ggplot(metadata, aes(x = log10(as.numeric(nFeature_RNA)+1)) )  +
       geom_histogram(aes(color = orig.ident, fill = orig.ident),  position = "identity", bins = 30, alpha = 0.4 )  + facet_grid(. ~ orig.ident) +
  xlab("") +
  ylab("log10(nFeature_RNA+1)") +
  ggtitle("Histogram ::  Genes(Features) per cell")

png(file = glue("{base.dir}/plots/{status}_{condition}_count_and_genes_per_cell_histogram.png"),width = 1500,height = 500)
print(counts.per.cell / genes.per.cell)
dev.off()


ordered.gene <- metadata[order(metadata$nFeature_RNA),] %>% ggplot(aes( x = seq_along(nFeature_RNA), y=nFeature_RNA, color = orig.ident)) + geom_point()  + facet_grid(. ~ orig.ident) +
  xlab("") +
  ylab("NFeatures") +
  ggtitle("Cells ordered by Nfeatures")

png(file = glue("{base.dir}/plots/{status}_{condition}_cells_ordered_by_Nfeatures.png"),width = 1500,height = 500)
print(ordered.gene)
dev.off()
    
###################################       Seurat QC        ############################################
    
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

density.plot.log10GenesPerCount <- metadata %>% 
  ggplot(aes(color=orig.ident, x = log10GenesPerCount, fill = orig.ident)) + 
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

  
density.plot.riboRatio <- metadata %>% 
  ggplot(aes(color=orig.ident, x = riboRatio, fill = orig.ident)) + 
   geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") 
    
png(file = glue("{base.dir}/plots/{status}_{condition}_RNA_density.png"),width = 1850,height = 500)
print(density.plot.nCount_RNA + density.plot.nFeature_RNA  + density.plot.log10GenesPerCount  + density.plot.mitoRatio + density.plot.riboRatio )
dev.off()
    
# -------------------------
png(file=glue("{base.dir}/plots/{status}_{condition}_RNA_vs_nFeature_RNA_scatter.png"),width = 1000,height = 500)
p <- FeatureScatter(object = data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
print(p)
dev.off()

png(file=glue("{base.dir}/plots/{status}_{condition}_nCount_RNA_vs_percent.mito_scatter.png"),width = 1000,height = 500)
p <- FeatureScatter(object = data, feature1 = "nCount_RNA", feature2 = "mitoRatio")
print(p)
dev.off()
    
png(file=glue("{base.dir}/plots/{status}_{condition}_nCount_RNA_vs_percent.ribo_scatter.png"),width = 1000,height = 500)
p <- FeatureScatter(object = data, feature1 = "nCount_RNA", feature2 = "riboRatio")
print(p)
dev.off()
    
    
print("::: QC PLOT - Seurat :::")

png(file=glue("{base.dir}/plots/{status}_{condition}_violonQC.png"),width = 1500,height = 500)
print(VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio","riboRatio","log10GenesPerCount"), ncol = 4)  + theme(  axis.text = element_text( size = 14)))
dev.off()

# Visualize the distribution of cell cycle markers across
png(file=glue("{base.dir}/plots/{status}_{condition}_cellcycle_ridgleplot.png"),width = 1500,height = 500)
print(RidgePlot(data, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2))
dev.off()

    
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
#                                  ident    sum    detected  total
# CTL_Vertes_AAACGAAAGACGCAGT-1 CTL_Vertes 22036     4927 22036

rowData(seurat.to.sce) <- cbind(rowData(seurat.to.sce),perFeatureQCMetrics(seurat.to.sce)) # mean / detected
# mean is the mean UMI/count per cell.
# detected is the number of umi/count for this feature.

print(head(colData(seurat.to.sce)))
print(head(rowData(seurat.to.sce)))

###################################      Scatter QC            #########################################

print ("::: QC PLOT - scater :::")

# Plot an overview of expression for each cell # Scater fonction...
png(file = glue("{base.dir}/plots/{status}_{condition}_n_cells_by_counts.png"),width = 1000,height = 500)
print(plotRowData(seurat.to.sce, y="detected", x="mean") +  scale_x_log10())
dev.off()

# Plot the percentage of counts accounted for by the top most highly expressed features across the dataset.
# Each row on the plot corresponds to a feature and is sorted by average expression
# These ticks can be coloured according to cell-level metadata, as specified by colour_cells_by

png(file=glue("{base.dir}/plots/{status}_{condition}_scater_mostexpresed.png"),width = 1600,height = 750)

p1 <- plotHighestExprs(seurat.to.sce[, seurat.to.sce$ident==condition], exprs_values = "counts",colour_cells_by="detected")   +  theme(  axis.text = element_text( size = 14), legend.text = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_blank(),legend.title = element_text( size = 14 )  ) + 
 xlab("Percentage of counts") + ylab("Genes") +labs(fill = "Counts per Cell")

print(p1)
      
dev.off()

return (seurat.to.sce)


 
}
                                              
#######################################################################################################
###################################       PCA            ########################################
#######################################################################################################

pca.stuffs <- function  (data,status,condition){
                                                 
# Perform PCA 

#VizPCA(object = pbmc, pcs.use = 1:2)

if (FALSE) { 
    data <- JackStraw(data, num.replicate = 100)
    data <- ScoreJackStraw(data, dims = 1:50)

    png(file = glue("{base.dir}/plots/{status}_{condition}_dimension_pca.png"),width = 1500,height = 500)
    j <- JackStrawPlot(data, dims = 1:50)
    e <- ElbowPlot(data.norm)
    print(e+j)
    dev.off() 
}
    
    print("Seurat Phase PCA")
    data <- RunPCA(data, features = VariableFeatures(data))

    png(file=glue("{base.dir}/plots/{status}_{condition}_genesBigVariance_pca.png"),width = 1500,height = 500)
    print(DimPlot(data,shape.by="orig.ident",group.by = "Phase",raster = FALSE,pt.size = 2))
    dev.off() 

    data <- RunPCA(data, features = c(cc.genes$s.genes, cc.genes$g2m.genes))

    png(file=glue("{base.dir}/plots/{status}_{condition}_genesCellCycle_pca.png"),width = 1500,height = 500)
    print(DimPlot(data,shape.by="orig.ident",group.by = "Phase",raster = FALSE,pt.size = 2) )   
    dev.off() 

    # ----------------------------------------------------------------------------------
    #Erreur : Cannot find 'phases.cyclone' in this Seurat object
    #print("Cyclone PCA")
    #data <- RunPCA(data, features = VariableFeatures(data))

    #png(file=glue("{base.dir}/plots/{status}_{analysed.conditions.string}_genesBigVariance_pca_cyclone.png"),width = 1500,height = 500)
    #print(DimPlot(data,shape.by="orig.ident",group.by = "phases.cyclone",raster = FALSE,pt.size = 2))
    #dev.off() 

    #data <- RunPCA(data, features = c(cc.genes$s.genes, cc.genes$g2m.genes))

    #png(file=glue("{base.dir}/plots/{status}_{analysed.conditions.string}_genesCellCycle_pca_cylone.png"),width = 1500,height = 500)
    #print(DimPlot(data,shape.by="orig.ident",group.by = "phases.cyclone",raster = FALSE,pt.size = 2) )   
    #dev.off() 
                          

#print("::: Data after Normalization & PCA :::")# data@meta.data
#print(head(data[[]]))# data@meta.data

}

#######################################################################################################
#######################################################################################################
###################################      MAIN             #############################################
#######################################################################################################
#######################################################################################################   
                                                     
                                                     
#######################################################################################################
###################################       Load h5 files        ########################################
#######################################################################################################
                                                     
plotting = opt$plot
set.seed(100)
# nCount_RNA = sum =  the number of UMIs per cell nUMI
# nFeature_RNA = detected ,the number of genes detected per cell nGene

# mean is the mean UMI/count per cell.
# detected is the number of umi/count for this feature

base.dir = opt$directory

#base.dir = "/data/villemin/data/toulouse/scRNAseqCells-1/CellRanger"
                                                     
dir.create(glue("{base.dir}/plots"), showWarnings = F)

experiments.list           <- as.list(args)
names(experiments.list)    <- as.list(experiments.list)
analysed.conditions.string <- paste(names(experiments.list), collapse = '.')

print(analysed.conditions.string)
#str(experiments.list)
                                                     
for (condition in experiments.list ) {
   
    
    # min.cells Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
    # min.features Include cells where at least this many features are detected.
    
    # Read 10x output & Create individual Seurat objects
experiments.list[[condition]]$seurat.object  <-  CreateSeuratObject(Read10X_h5(glue("{base.dir}/{condition}/outs/filtered_feature_bc_matrix.h5")), project = condition, assay = "RNA", min.cells = 3, min.features = 200) 
  

   experiments.list[[condition]]$seurat.object$log10GenesPerCount <-  log10(experiments.list[[condition]]$seurat.object$nFeature_RNA) / log10(experiments.list[[condition]]$seurat.object$nCount_RNA)
   # Compute percent mito ratio & riboRatio

   experiments.list[[condition]]$seurat.object$mitoRatio          <-  PercentageFeatureSet(object = experiments.list[[condition]]$seurat.object, pattern = "^MT-")
   experiments.list[[condition]]$seurat.object$riboRatio          <-  PercentageFeatureSet(object = experiments.list[[condition]]$seurat.object, pattern = "^RP[SL]")
   experiments.list[[condition]]$dimension                        <-  dim(x = experiments.list[[condition]]$seurat.object)
    
    print("====> Control.Quality")
    # Create metadata dataframe
    seurat.to.sce <- control.quality (experiments.list[[condition]]$seurat.object,"Raw",condition)
    
    #map2( ExamResults, MathsRemark, assign_in, where=list("Maths",2) )

    print("====> Cell.cycle")
    # Normalisation is  done here
    experiments.list[[condition]]$seurat.object <- cell.cycle (experiments.list[[condition]]$seurat.object,"Raw",condition,seurat.to.sce)
    
   
}
                                                     
###############################################################################################################################################
###################################        Merge        #######################################################################################
###############################################################################################################################################
             
if (length(experiments.list) > 1 ) {  
    
print("::: Merge :::")
#By default, merge() will combine the Seurat objects based on the raw count matrices, erasing any previously normalized and scaled data matrices
#   merge.data = TRUE	
# Merge the data slots instead of just merging the counts (which requires renormalization);
# this is recommended if the same normalization approach was applied to all objects

data <- merge(x = experiments.list[[1]]$seurat.object,y = unlist(map(experiments.list[-1],`[`,c("seurat.object")), use.names=FALSE), merge.data = TRUE,add.cell.ids = names(experiments.list))

table(data$orig.ident)

head(data@meta.data)
                                                     
###############################################################################################################################################
###################################        FindVariableFeatures & Scaled (Normalize has been done previously)     #############################       
###############################################################################################################################################
                                                     
#print ("::: Normalized merge data :::")
#GetAssayData(pbmc.combined)[1:10, 1:15]
#data <- NormalizeData(data) if no merge.data = TRUE in merge
                                                     
print ("::: FindVariableFeatures & Scale :::")
                    
# These part need to be applied on a merge object.

#data <- NormalizeData(data) No need to normalise if  merge.data = TRUE is set in merge
#data <- FindVariableFeatures(data, selection.method = "vst",nFeatures = 2000) #nFeatures=2000
#data <- ScaleData(data, features = rownames(data)) # After normalisation

print ("::: PCA stuffs :::")
    
#pca.stuffs(data,"Raw",analysed.conditions.string)
                         
saveRDS(data,file = glue("{base.dir}/{analysed.conditions.string}_raw.rds"))
              
print ("::: End :::")

}