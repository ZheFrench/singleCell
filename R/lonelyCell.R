library(Seurat)
library(glue)
library(hdf5r)
library(scater) # Davis McCarthy Bioinformatics 2017 Scater: pre-processing, quality control, normalisation and visualisation of single-cell RNA-seq data in R
library(patchwork)
library(SingleCellExperiment)# Robert A. Amzequita 2019 Nature merthods Orchestrating single-cell analysis with Bioconductor
library(scran)
library(dplyr)
library(org.Hs.eg.db)

#######################################################################################################
###################################       Load h5 files        ########################################
#######################################################################################################

set.seed(100)
# nCount_RNA = the number of UMIs per cell nUMI
# nFeature_RNA = the number of genes detected per cell nGene

cond1 = "OSI_TIPI_Vertes"
cond2 = "OSI_TIPI_Rouges"

base.dir = "/data/villemin/data/toulouse/scRNAseqCells-2/CellRanger"
dir.create(glue("{base.dir}/plots"))

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

#metadata %>% ggplot(aes( x = ?, color = Cell_line,fill = group)) + geom_density(alpha = 0.2, fill = "green") + scale_x_log10() + theme_classic() )
data <- merge(cond1.object, cond2.object, add.cell.ids=c("OSI_TIPI_Vertes","OSI_TIPI_Rouges"))
#head(data@meta.data)
View(data@meta.data)

# -------------------------
# Create metadata dataframe
metadata <- data@meta.data

metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^OSI_TIPI_Vertes_"))] <- "OSI_TIPI_Vertes"
metadata$sample[which(str_detect(metadata$cells, "^OSI_TIPI_Rouges_"))] <- "OSI_TIPI_Rouge"
data@meta.data <- metadata

View(data@meta.data)

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_nCount_RNA_density.png"),width = 500,height = 500)

# Visualize the number UMIs/transcripts per cell
metadata %>% 
  ggplot(aes(color=sample, x=nCount_RNA, fill= sample)) + 
   geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)
dev.off()
dev.off()
stop()
# -------------------------


#######################################################################################################
###################################       Seurat QC        ############################################
#######################################################################################################
print ("::: QC PLOT - Seurat :::")

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_nFeature_RNA.png"),width = 500,height = 500)
VlnPlot(data, features = "nFeature_RNA", pt.size = 0.1) + NoLegend() + theme(  axis.text = element_text( size = 14))
dev.off()

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_nCount_RNA.png"),width = 500,height = 500)
VlnPlot(data, features = "nCount_RNA", pt.size = 0.1) + NoLegend() + theme(  axis.text = element_text( size = 14))
dev.off()

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_percent.Ratio.png"),width = 500,height = 500)
VlnPlot(data, features = "mitoRatio", pt.size = 0.1) + NoLegend() + theme(  axis.text = element_text( size = 14))
dev.off()

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_log10GenesPerUMI.png"),width = 500,height = 500)
VlnPlot(data, features = "log10GenesPerUMI", pt.size = 0.1) + NoLegend() + theme(  axis.text = element_text( size = 14))
dev.off()

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_violonQC.png"),width = 1500,height = 500)
VlnPlot(object = data, features = c("nFeature_RNA", "nCount_RNA", "mitoRatio","log10GenesPerUMI"), ncol = 4)  + theme(  axis.text = element_text( size = 14))
dev.off()

#######################################################################################################
###################################       SingleCellExperiment QC   (scater)     ######################
#######################################################################################################

# Switch seurat object to SingleCellExperiment (Nat methods 2020)
seurat.to.sce <- as.SingleCellExperiment(data)

# Check the assays present
#assays(seurat.to.sce)
#head(seurat.to.sce)
if (FALSE) { 
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

print ("::: QC PLOT - scater :::")

# Plot an overview of expression for each cell
png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_n_cells_by_counts.png"),width = 1000,height = 500)
plotRowData(seurat.to.sce, y="detected", x="mean") +  scale_x_log10()
dev.off()

# Plot the percentage of counts accounted for by the top most highly expressed features across the dataset.
# Each row on the plot corresponds to a feature and is sorted by average expression
# These ticks can be coloured according to cell-level metadata, as specified by colour_cells_by

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_scater_mostexpresed.png"),width = 1600,height = 750)

p1 <- plotHighestExprs(seurat.to.sce[, seurat.to.sce$ident=="OSI_TIPI_Rouges"], exprs_values = "counts",colour_cells_by="detected")   +  theme(  axis.text = element_text( size = 14), legend.text = element_text(size=14), axis.title.x = element_text(size=14), axis.title.y = element_blank(),legend.title = element_text( size = 14 )  ) + 
 xlab("Percentage of counts") + ylab("Genes") +labs(fill = "Counts per Cell")

p2 <- plotHighestExprs(seurat.to.sce[, seurat.to.sce$ident=="OSI_TIPI_Vertes"], exprs_values = "counts",colour_cells_by="detected")  + theme( axis.text = element_text( size = 14), legend.text = element_text(size=14), axis.title.x = element_text(size=14),axis.title.y = element_blank(),legend.title = element_text( size = 14 ) ) +
xlab("Percentage of counts") + ylab("Genes") + guides(fill=guide_legend(title="Counts per Cell"))
print(p1 + p2 + plot_layout(ncol=2) )

dev.off()
}#theme(base_size = 22)
#######################################################################################################
###################################       Calculate cell-cycle scores            ######################
#######################################################################################################

hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))



organism = "hsapiens_gene_ensembl"
host="mar2016.archive.ensembl.org" # ensembl 84 ... There is no 83 ...
symbol_description='hgnc_symbol'


# BIOMART OBJECT
edb = useMart("ENSEMBL_MART_ENSEMBL", dataset=organism,host=host)
colnames(sampleTable)[colnames(sampleTable) == 'gene'] <- 'ensembl_gene_id'
colnames(sampleTable)[colnames(sampleTable) == 'Gene'] <- 'ensembl_gene_id'



# Retrieve gene infos And entrezeneId needed fr KEGGPATHWAY
gene_infos = getBM(attributes=c('ensembl_gene_id',symbol_description,'gene_biotype'),values=sampleTable$ensembl_gene_id,filters='ensembl_gene_id',mart=edb)
# Correction of bug 
# Output had duplicate lines
gene_infos_without_dup <- gene_infos[!duplicated(gene_infos$ensembl_gene_id),]


res_annotated <- join(gene_infos, sampleTable, by='ensembl_gene_id', type='left', match='all')

anno <- select(org.Hs.eg.db, keys=rownames(seurat.to.sce),keytype="SYMBOL", column="ENSEMBL")

ensembl <- anno$ENSEMBL[ match(rownames(seurat.to.sce), anno$SYMBOL)]
# Use only genes related to biological process cell cycle to speed up
# https://www.ebi.ac.uk/QuickGO/term/GO:0007049 = cell cycle (BP,Biological
# Process)
GOs <- na.omit(select(org.Hs.eg.db, keys = na.omit(ensembl), keytype = "ENSEMBL", column = "GO"))
## 'select()' returned many:many mapping between keys and columns
GOs <- GOs[GOs$GO == "GO:0007049", "ENSEMBL"]

GOs <- GOs[GOs$GO == "GO:0007049", "ENSEMBL"]
hs.pairs <- lapply(hs.pairs, function(x) {
    x[rowSums(apply(x, 2, function(i) i %in% GOs)) >= 1, ]
})
str(hs.pairs)
cc.ensembl <- ensembl[ensembl %in% GOs]

# Using Ensembl IDs to match up with the annotation in 'hs.pairs'.
assignments <- cyclone(seurat.to.sce[ensembl %in% cc.ensembl, ], hs.pairs, ensembl[ensembl %in% cc.ensembl])

png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_cellcyle.png"),width = 1500,height = 500)
plot(assignments$score$G1, assignments$score$G2M,xlab="G1 score", ylab="G2/M score", pch=16)
dev.off()

#seurat.to.sce$G1.score <- assignments$scores$G1
#seurat.to.sce$G2M.score <- assignments$scores$G2M
#seurat.to.sce$S.score <- assignments$scores$S
              
table(assignments$phases, (seurat.to.sce))
      
png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_n_cells_by_counts.png"),width = 1000,height = 500)
      
p1 <- plotColData(seurat.to.sce, y = "G2M.score", x = "G1.score", colour_by = "sample")
p2 <-plotColData(seurat.to.sce, y = "G2M.score", x = "sample", colour_by = "sample")
p3 <-plotColData(seurat.to.sce,y = "G1.score", x = "sample", colour_by = "sample")
p4 <- plotColData(seurat.to.sce, y = "S.score", x = "sample", colour_by = "sample")
p1 + p2 + p3 + p4 
              
dev.off()      

              
seurat.to.sce <- runPCA(seurat.to.sce, ncomponents=20)
png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_pca.png"),width = 1000,height = 500)
plotPCA(seurat.to.sce, ncomponents = 4, colour_by = "level1class")              
dev.off()
              
stop("kikou")
              
png(file=glue("{base.dir}/plots/{cond1}_vs_{cond2}_PCA.png"),width = 1500,height = 500)
plotPCA(seurat.to.sce, colour_by=I(assignments$phases), shape_by=I(relabel)) +  ggtitle("Before") 
dev.off()
      
      


#Low library size: When cells are very degraded or absent from the library preparation, the number of reads sequenced from that library will be very low. It’s important to remove these cells from downstream analyses.

#Low number of expressed genes: A low number of expressed genes may be a result of poor-quality cells (e.g. dying, degraded, damaged, etc.), followed by high PCR amplification of the remaining RNA. Again, these cells should be removed from downstream analyses.

#High mitochondrial gene content: High concentrations of mitochondrial genes is often a result of damaged cells where the endogenous RNA escapes or degrades. As mitochondria has its own cell membranes, it is often the last DNA/RNA in damaged cells to degrade and hence occurs in high quantities during sequencing.

#Batch effect: Large scRNA-seq projects usually need to generate data across multiple batches due to logistical constraints. However, the processing of different batches is often subject to variation, e.g., changes in operator, differences in reagent quality and concentration, the sequencing machine used, etc. This results in systematic differences in the observed expression in cells from different batches, which we refer to as “batch effects”. Batch effects are problematic as they can be major drivers of variation in the data, masking the relevant biological differences and complicating interpretation of the results.
# Filter out low quality reads using selected thresholds - these will change with experiment

# 4006_red
filtered_seurat.cond1 <- subset(x = cond1.object,subset= (nUMI >= 9000) & (nGene >= 2800) & (log10GenesPerUMI >= 0.8) & (mitoRatio < 0.20))
# 4006_green
filtered_seurat.cond2 <- subset(x = cond2.object, subset= (nUMI >= 9000) & (nGene >= 2800) & (log10GenesPerUMI >= 0.8) & (mitoRatio < 0.20))

seurat <- CreateSeuratObject(counts = cbind(data_condition1, data_condition2), min.cells = 3, min.features  = 200, project = cbind("OSI_TIPI_Vertes", "OSI_TIPI_Rouges"), assay = "RNA")

head(seurat@meta.data)

## Explore the raw counts for the dataset
#dim(counts(sce))

#counts(sce)[1:6, 1:6]

#pbmc@meta.data$stim <- c(rep("STIM", ncol(stim.sparse)), rep("CTRL", ncol(ctrl.sparse)))

#p1 <- plotExpression(pbmc.sce, features = "MS4A1", x = "ident") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
#p2 <- plotPCA(pbmc.sce, colour_by = "ident")
#p1 + p2
#filtered_seurat.cond1 <- subset(x = cond1.object,subset= (nUMI >= 9000) &(nGene >= 2800) &(log10GenesPerUMI >= 0.8) & (mitoRatio < 0.20))

stop()