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
# speudoBulk.R 
#
# Source : 
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/OSCABioc2019__OSCABioc2019/
# http://biocworkshops2019.bioconductor.org.s3-website-us-east-1.amazonaws.com/page/muscWorkshop__vignette/
# https://hbctraining.github.io/scRNA-seq/lessons/pseudobulk_DESeq2_scrnaseq.html
# https://satijalab.org/seurat/archive/v3.0/de_vignette.html
# https://github.com/csoneson/conquer_comparison/tree/master/scripts
# https://f1000research.com/articles/5-2122/v2 next step when you want to compare one vs all others.
#
# Description : 
#
# SpeudoBulk.R differential Epresssion (DE) of scRNA-SEQ
# Randomly assign equal number of cells in 3 groups, sum up reads, and do differential expression analysis
#
#################################################################

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
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(Matrix.utils))
suppressPackageStartupMessages(library(BiocParallel))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(pheatmap))
suppressPackageStartupMessages(library(RColorBrewer))

#######################################################################################################
###################################       Load h5 files        ########################################
#######################################################################################################
plotting = FALSE
set.seed(100)
# nCount_RNA   = the number of UMIs per cell nUMI
# nFeature_RNA = the number of genes detected per cell nGene
# ------------------------------------------------------------------------------------------------------------
#cond1 = "OSI_TIPI_Vertes"
#cond2 = "OSI_TIPI_Rouges"
#cond1 = PC9_verte
#cond2 = PC9_rouge
#cond1 = "4006_verte"
#cond2 = "4006_rouge"
#base.dir = "/data/villemin/data/toulouse/scRNAseqCells-1/CellRanger"
cond1 = "CTL_Vertes"
cond2 = "CTL_Rouges"

base.dir = "/data/villemin/data/toulouse/scRNAseqCells-2/CellRanger"

# ------------------------------------------------------------------------------------------------------------
dir.create(glue("{base.dir}/DE"), showWarnings = F)

    
print(" Loading RDS.. ")
    
seurat.Object <- readRDS(file = glue("{base.dir}/{cond1}_vs_{cond2}.rds"))
print(head(seurat.Object[[]])) # data@meta.data

################################################################################################################
###################################      Create Speudo-Bulk data          ########################################
################################################################################################################
 
# Switch seurat object to SingleCellExperiment (Nat methods 2020)
sce.Object       <- as.SingleCellExperiment(seurat.Object)
sce.Object$ident <- sce.Object$orig.ident # Important due to a bug.

print(glue("# Total Genes {dim(sce.Object)[1]}"))
# Remove lowly expressed gene at cell level
# Sum of each gene equal at least one UMI, and this criteria should be met amongst 10 cells. (per conditions)
sce.Object <- sce.Object[rowSums(counts(sce.Object) > 1) >= 10, ]
print(glue("# Filtered Genes {dim(sce.Object)[1]}"))

length(sce.Object$ident)#5775

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


#######################################################################################################
###################################      DIFF EDGER          ##########################################
#######################################################################################################

print("##### INITIAL SINGLE CELLS ######")

#print("coldData")
#print(head(colData(sce.Object)))

#print("rowData")
#print(head(rowData(sce.Object)))

sce.Object$sample_id    <- paste(sce.Object$orig.ident, "-", sce.Object$subgroup,sep = "")

# Not clear what kind of  witchcraft is used with m & subset.exp . Finger crossed.
m          <- match(unique(sce.Object$sample_id), sce.Object$sample_id)
n_cells    <- as.numeric(table(sce.Object$sample_id))
subset.exp <- data.frame(colData(sce.Object)[m, ], n_cells, row.names = NULL) 
n_cells.speudoBulk.resume   <- subset.exp[, c("sample_id","orig.ident","n_cells")]
print(n_cells.speudoBulk.resume)
write.table(n_cells.speudoBulk.resume ,file = glue("{base.dir}/DE/{cond1}_{cond2}-nCell-speudoBulks.tsv"),quote=F,row.names=F,sep="\t")

# Aggregate across cluster-sample groups
groups <- colData(sce.Object)[, c("sample_id")]
# split by cluster, transform & rename columns
matrix.rnaseq <- aggregate.Matrix(t(counts(sce.Object)), groupings = groups, fun = "sum")   # Sum can be biased if you have a big diff in number of cells between two groups  

counts=t(matrix.rnaseq)
# ---------------------------------------------------------------------------------------------
# Remove lowly expressed gene at speudo-bulk level
# At leat 5 reads in 3 of pseudo bulk samples over 6

good <- apply(counts,1,function(x) sum(x>5))>=3
print(glue("# Filtered Genes {sum(good)}"))
counts <- counts[good,]

# Remove  method from JC based on TMM, no need to do that twice, edgeR take into account raw reads.
# https://www.biostars.org/p/317701/ 
# q <- apply(counts,2,function(x) quantile(x[x>0],prob=0.75))
# ncounts <- sweep(counts,2,q/median(q),"/")

# ---------------------------------------------------------------------------------------------

           
# construct SCE of pseudo-bulk counts
matrix.rnaseq.cse <- SingleCellExperiment(assays =list(counts=counts))
# Add QC metrics from scater
colData(matrix.rnaseq.cse) <- cbind(colData(matrix.rnaseq.cse),perCellQCMetrics(matrix.rnaseq.cse)) # sum / detected / total , threshold > 0 at least one count
rowData(matrix.rnaseq.cse) <- cbind(rowData(matrix.rnaseq.cse),perFeatureQCMetrics(matrix.rnaseq.cse)) # mean / detected is %

print("##### SPEUDO BULK ######")

print("coldData")
print(head(colData(matrix.rnaseq.cse)))
#write.csv(colData(matrix.rnaseq.cse), file=glue("{base.dir}/DE/{cond1}_{cond2}_colData.csv"), row.names = TRUE)

print("rowData")
print(head(rowData(matrix.rnaseq.cse)))
#write.csv(rowData(matrix.rnaseq.cse), file=glue("{base.dir}/DE/{cond1}_{cond2}_rowData.csv"),  row.names = TRUE)


#######################################################################################################
###################################      DIFF EDGER          ##########################################

#######################################################################################################
sce <- logNormCounts(matrix.rnaseq.cse) 
#var <- modelGeneVar(sce, BPPARAM = MulticoreParam(workers = 8), assay.type = "logcounts")

dge <- convertTo(sce, type="edgeR") #DESeq2

dge$samples$group <- sub("^(.*)-[0-9]", "\\1", rownames(dge$samples))
dge$samples$group <- relevel(factor(dge$samples$group),ref=cond2)

comp     <- glue("{cond1}-{cond2}")

de.design <- model.matrix(~0 + dge$samples$group)
colnames(de.design) <- gsub("^dge\\$samples\\$group","",colnames(de.design))

cm     <- makeContrasts(contrasts = comp,levels = dge$samples$group)

dge    <- estimateDisp(dge, de.design,robust=T)

fit.y  <- glmFit(dge, de.design)
lrt    <- glmLRT(fit.y,contrast=cm)

head(de.design)
#names(dge)
#names(lrt)
#######################################################################################################
###################################      FILES          ###############################################
#######################################################################################################
    
print("Writing differential results...")

result <- as.data.frame( topTags(lrt, adjust.method="BH",n=Inf, sort.by="PValue", p.value=1))
result <- cbind(genes = rownames(result), result, row.names = NULL)

result     <- cbind(result,dge$counts[result$genes,])


write.table(result,file = glue("{base.dir}/DE/{cond1}_{cond2}-differential.tsv"),quote=F,row.names=F,sep="\t")

result_up    <- subset(result, ( logFC >= 1.5 & FDR < 0.05 ))
result_up <- result_up[order(abs(result_up$logFC),decreasing = TRUE),]
write.table(result_up,file = glue("{base.dir}/DE/{cond1}_{cond2}-differential-up.tsv"),quote=F,row.names=F,sep="\t")

result_down <- subset(result, ( logFC <= -1.5 & FDR < 0.05 ))
result_down <- result_down[order(abs(result_down$logFC),decreasing = TRUE),]
write.table(result_down,file = glue("{base.dir}/DE/{cond1}_{cond2}-differential-down.tsv"),quote=F,row.names=F,sep="\t")

summary <- as.data.frame(summary(dt<-decideTestsDGE(lrt, adjust.method="BH",p.value = 0.05,lfc = 1.5)))
colnames(summary)[1] <- "FC"
colnames(summary)[2] <- "Analyse"
write.table(summary ,file = glue("{base.dir}/DE/{cond1}_{cond2}-summary.tsv"),quote=F,row.names=F,sep="\t")

#######################################################################################################
###################################      PLOTS          ###############################################
#######################################################################################################
    
print("Plotting...")
#Plot the genewise biological coefficient of variation (BCV) against gene abundance (in log2 counts per million).
png(file = glue("{base.dir}/DE/{cond1}_{cond2}_BCV.png"),width = 500,height = 500)
plotBCV(dge,xlab="Average log CPM", ylab="Biological coefficient of variation", pch=16, cex=0.2) 
dev.off()

#Both of these functions plot the log-fold change (i.e. the log of the ratio of expression levels for each gene between two experimential groups) against the log-concentration (i.e. the overall average expression level for each gene across the two groups). To represent counts that were low (e.g. zero in 1 library and non-zero in the other) in one of the two conditions, a 'smear' of points at low A value is presented in 'plotSmear'
de.genes <- rownames(topTags(lrt, adjust.method="BH",n=Inf, sort.by="PValue", p.value=0.05)$table)
png(file = glue("{base.dir}/DE/{cond1}_{cond2}_smear.png"),width = 500,height = 500)
plotSmear(lrt, de.tags=de.genes,xlab="Average logCPM", ylab="logFC", pch=19, cex=0.2)
dev.off()

#Plot samples on a two-dimensional scatterplot so that distances on the plot approximate the typical log2 fold changes between the samples.
png(file = glue("{base.dir}/DE/{cond1}_{cond2}_mds.png"),width = 500,height = 500)
limma::plotMDS(dge, col = as.numeric(as.factor(dge$samples$group)) , labels= dge$samples$group , pch = 24)
dev.off()

## Obtain logical vector regarding whether padj values are less than 0.05
threshold_OE <- result$FDR < 0.05 
## Determine the number of TRUE values
length(which(threshold_OE))
## Add logical vector as a column (threshold) to the res_tableOE
result$threshold <- threshold_OE 

## Sort by ordered logFC
result <- result[order(abs(result$logFC), decreasing = TRUE), ] 
## Create a column to indicate which genes to label
n.top = 20
result$genelabels <- ""
result$genelabels[1:n.top] <- result$genes[1:n.top]
IDS <- as.character(result$genelabels[1:n.top])

## Volcano plot log10(0.1)-> -1
png(file = glue("{base.dir}/DE/{cond1}_{cond2}_volcano.png"),width = 500,height = 1000)
ggplot(as.data.frame(result)) +
        geom_point(aes(x = logFC, y = -log10(FDR + 0.1), colour = threshold)) +
        geom_text_repel(aes(x = logFC, y = -log10(FDR + 0.1), label = ifelse(as.character(result$genes) %in% IDS  , IDS,"")),max.overlaps = Inf) +
        xlab("Log2 fold change") + 
        ylab("-Log10 FDR+0.1")        +
        geom_hline(yintercept=-log10(0.1+0.05) ) + 
        geom_vline(xintercept=c(-1.5,1.5) , linetype=c("dashed","dashed"))  +
        theme(legend.position = "none" ) 
dev.off()

### Run pheatmap using n.top
df.heatmap <- result[result$genes %in% IDS,]
rownames(df.heatmap) <- df.heatmap$genes
drop <- c("genelabels","threshold")
df.heatmap <-df.heatmap[,!(names(df.heatmap) %in% drop)]   
### Annotate our heatmap (optional)
#annotation <- data.frame(sampletype=mov10_meta[,'sampletype'],row.names=rownames(mov10_meta))
#annotation = annotation,anno <- as.data.frame(colData(vsd)[, c("cell","dex")])
### Set a color palette
heat_colors <- brewer.pal(6, "YlOrRd")
png(file = glue("{base.dir}/DE/{cond1}_{cond2}_heatmap.png"),width = 500,height = 1000)
pheatmap(df.heatmap[,9:length(names(df.heatmap)) ], color = heat_colors, cluster_rows = T, show_rownames=T,border_color=NA, fontsize = 10, scale="row",fontsize_row = 10, height=20)
dev.off()

