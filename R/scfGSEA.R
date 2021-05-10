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
  make_option(c("-d", "--directory"), type="character", default=NULL, help="Directory to look for output directories of CellRanger. ", metavar="PATH2DIRECTORY")

)

parser    = OptionParser(usage = "%prog [options] file ",option_list=option_list);
arguments = parse_args(parser, positional_arguments = TRUE);
opt  <- arguments$options
args <- arguments$args

print("> OPTS : ")
print(opt$directory)
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
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tibble))

suppressPackageStartupMessages(library(fgsea))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(reshape2))


#http://barc.wi.mit.edu/education/hot_topics/scRNAseq_2020/SingleCell_Seurat_2020.html


#######################################################################################################
#######################################################################################################
###################################      MAIN             #############################################
#######################################################################################################
#######################################################################################################   



base.dir = opt$directory
    

dir.create(glue("{base.dir}/GSEA"), showWarnings = F)

#https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
#https://github.com/ctlab/fgsea
#https://stephenturner.github.io/deseq-to-fgsea/



set.seed(423)
#------------------------------------------------------------------------------
#   Read directory with files previously created toorder by FC and apply fsGEA
#------------------------------------------------------------------------------
files.list <- list.files(glue("{base.dir}"),pattern="(*)-differential.tsv$")

#h.All <- gmtPathways("/data/villemin/annotation/gsea/MSigDB/h.all.v7.2.symbols.gmt")# 50
#h.All.bis <- read.gmt("/data/villemin/annotation/gsea/MSigDB/h.all.v7.2.symbols.gmt")
file.gmt <- "/data/villemin/annotation/gsea/MSigDB/h.all.v7.2.symbols.gmt"
#file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c5.go.v7.2.symbols.gmt"

#file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c2.all.v7.2.symbols.gmt"
#file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c3.tft.gtrd.v7.2.symbols.gmt"
#file.gmt <- "/data/villemin/annotation/gsea/MSigDB/c6.all.v7.2.symbols.gmt"
#file.gmt <- "/data/villemin/annotation/gsea/MSigDB/jp.gmt"
database <- basename(file.gmt)

final.file.padj <- glue("{base.dir}/GSEA/fgsea.{database}.padj.txt")
final.file.nes  <- glue("{base.dir}/GSEA/fgsea.{database}.nes.txt")


h.All.bis <- gmtPathways(file.gmt) # 6226
h.All <- read.gmt(file.gmt)

subDir = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file.gmt))

i = 0
colname <-"pathway"
for (file in files.list){
       
  asbolutepath2file <- glue("{base.dir}/{file}")

  file <- gsub(file,pattern= "-differential.tsv",replacement="",perl=T)

  subpath <- file.path(glue("{base.dir}/GSEA/{file}/"), subDir)
  dir.create(subpath, recursive = TRUE, showWarnings = F)
  
  print(file)
  print(subpath)

  dataframe.expression <- fread(asbolutepath2file,data.table=F)

  dataframe.expression <- subset(dataframe.expression,select=c(genes,logFC))

  dataframe.expression <- dataframe.expression[order(dataframe.expression$logFC),]

  # ClusterProfiler need a decreasing order...
  dataframe.expression.decreasing <- dataframe.expression[order(dataframe.expression$logFC, decreasing = TRUE),]
  ranks_decreasing <- deframe(dataframe.expression.decreasing)
  
  ranks <- deframe(dataframe.expression)

  #not using fgseaMultilevel...a really tiny difference in NES score due to the fact it use fgsea.
  # But genes used are the same so I used it to plot with heatplot function of clusterProfiler the leading edge genes contained in the signature I am interested in
  egmt2 <- GSEA(ranks_decreasing, TERM2GENE = h.All, by = "fgsea",verbose = TRUE ,nPermSimple = 10000 ,minGSSize  = 10, maxGSSize  = 325 , eps = 0,  pvalueCutoff = 1)
  # You dont need the whole object to be written
  write.table(egmt2, file=glue("{subpath}/{file}-full-gsea-clusterprofiler.tsv"),quote=F,row.names=F,sep="\t")

  # Yeah I do it again (I know.Don't say a fucking word moron.)
  fgseaRes     <- fgseaMultilevel(pathways= h.All.bis, stats=ranks,eps = 0, nPermSimple = 10000 ,minSize  = 10, maxSize  = 325)
  
  for (pathway in names(h.All)){
    if (pathway %in% c("BLUM_RESPONSE_TO_SALIRASIB_DN","BLUM_RESPONSE_TO_SALIRASIB_UP" )){
      #print(pathway)
      #print(h.All[[pathway]])
      p<- plotEnrichment(h.All.bis[[pathway]], ranks) + labs(title=pathway)
    
      png(file=glue("{subpath}/{file}-{pathway}.png"))
      print(p)
      dev.off()
      
    }  
  }
  
  fgseaResTidy <- fgseaRes %>% as_tibble()

  fwrite( filter(fgseaResTidy ,padj <= 1), file=glue("{subpath}/{file}-full-fgsea.tsv"), sep="\t", sep2=c("", "/", ""))

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n = 10), pathway]

  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n = 10), pathway]

  
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  png(file=glue("{subpath}/{file}-top-global.png"),width = 900)
  plotGseaTable(h.All.bis[topPathways], ranks, fgseaRes, gseaParam = 0.5) 
  dev.off()
  
  
  png(file=glue("{subpath}/{file}-global.png"),width = 900)
  print(ggplot( filter(fgseaResTidy ,abs(NES) >= 1.7), aes(reorder(pathway, NES), NES)) +
    geom_col( aes(fill=padj<0.01),color="#5F6368") +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",title=file) + 
    theme_minimal())
  dev.off()
  
  test <- c("PRC2_EZH2_UP.V1_UP","PRC2_EZH2_UP.V1_DN","BLUM_RESPONSE_TO_SALIRASIB_DN","BLUM_RESPONSE_TO_SALIRASIB_UP","E2F1_UP.V1_UP","E2F1_UP.V1_DN")
  
  #write.table(dataframe.expression,glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}/{file}-expresssion.tsv"),quote=F,row.names=T,sep="\t")

  heatplot_up <- heatplot(egmt2,foldChange = ranks, showCategory =   topPathwaysUp)
  #png(file=glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}/{file}-heatplot-up-gsea.png"),width=2000)   
  #print(  heatplot_up + theme(  axis.text.y = element_text( size = 12)))
  #dev.off()

  heatmap.reshaped <- dcast(heatplot_up$data, categoryID~Gene, value.var='foldChange')
  
  rownames(heatmap.reshaped)<- heatmap.reshaped[,1]
  heatmap.reshaped<- heatmap.reshaped[,-1]
  heatmap.reshaped[is.na(heatmap.reshaped)] <- 100
  #heatmap.reshaped[heatmap.reshaped !=0] <- 1
  #write.table(heatplot_up$data,glue("{subpath}/expresssion-up-gsea.tsv"),quote=F,row.names=F,sep="\t")
  write.table(heatmap.reshaped,glue("{subpath}/{file}-genes-nes-up.tsv"),quote=F,row.names=T,col.names=NA,sep="\t")

  heatplot_down <-  heatplot(egmt2,foldChange=ranks, showCategory = topPathwaysDown)
  #png(file=glue("{subpath}/{file}-heatplot-down-gsea.png"),width=2000)   
  #print(heatplot_down + theme(  axis.text.y = element_text( size = 12))) 
  #dev.off()
  
  heatmap.reshaped <- dcast(heatplot_down$data, categoryID~Gene, value.var='foldChange')
  rownames(heatmap.reshaped)<- heatmap.reshaped[,1]
  heatmap.reshaped<- heatmap.reshaped[,-1]
  
  heatmap.reshaped[is.na(heatmap.reshaped)] <- 100
  #heatmap.reshaped[heatmap.reshaped !=0] <- 1
  
  #write.table(heatplot_down$data,glue("{subpath}/expresssion-down-gsea.tsv"),quote=F,row.names=F,sep="\t")

  write.table(heatmap.reshaped,glue("{subpath}/{file}-genes-nes-down.tsv"),quote=F,row.names=T,col.names=NA,sep="\t")
  
  final.padj <-  fgseaResTidy %>% select(-leadingEdge, -pval,-log2err,-size, -ES,-NES)%>% rename(!!file  :=  padj)
  final.nes  <-  fgseaResTidy %>% select(-leadingEdge, -pval,-log2err,-size, -ES,-padj)%>% rename(!!file :=  NES)
  
  colname[i] <- file
  
  if (!exists("dataset.padj") ){
    dataset.padj <- final.padj
    dataset.nes  <- final.nes
    next
    }
  
  # if the merged dataset does exist, append to it 
  if (exists("dataset.padj")){
    
    temp_dataset.padj <-  fgseaResTidy %>% select(-leadingEdge, -pval,-log2err,-size, -ES,-NES)%>% rename( !!file := padj)
    temp_dataset.nes  <-  fgseaResTidy %>% select(-leadingEdge, -pval,-log2err,-size, -ES,-padj)%>% rename(!!file :=  NES)
    
    dataset.padj <- inner_join(dataset.padj,temp_dataset.padj, by =c("pathway"),keep=FALSE)
    dataset.nes<- inner_join(dataset.nes,temp_dataset.nes, by =c("pathway"),keep=FALSE)
    
    rm(temp_dataset.padj)
    rm(temp_dataset.nes)
    
  }
i=i+1
}



dataset.padj <-  dataset.padj %>% select(-contains(".y"))
dataset.nes  <-  dataset.nes %>% select(-contains(".y"))

write.table(dataset.padj,file=final.file.padj,quote=F,row.names=F,sep="\t")
write.table(dataset.nes,file=final.file.nes,quote=F,row.names=F,sep="\t")