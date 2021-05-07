library(optparse)
library(fgsea)
library(data.table)
library(tidyr)
library(dplyr)
library(glue)
library(tibble)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(reshape2)
#https://davetang.org/muse/2018/01/10/using-fast-preranked-gene-set-enrichment-analysis-fgsea-package/
#https://github.com/ctlab/fgsea
#https://stephenturner.github.io/deseq-to-fgsea/


cellline = "pc9"

final.file.padj <- glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}-fgsea.padj.txt")
final.file.nes  <- glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}-fgsea.nes.txt")

set.seed(423)
#------------------------------------------------------------------------------
#   Read directory with files previously created toorder by FC and apply fsGEA
#------------------------------------------------------------------------------
files.list <- list.files(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}/"),pattern="(*)-fgsea-gene-ratios.txt$")

#h.All <- gmtPathways("/home/jp/eclipse-workspace/database/MSigDB/h.all.v7.2.symbols.gmt")# 50
#h.All.bis <- read.gmt("/home/jp/eclipse-workspace/database/MSigDB/h.all.v7.2.symbols.gmt")

#file.gmt <- "/home/jp/eclipse-workspace/database/MSigDB/c2.all.v7.2.symbols.gmt"
file.gmt <- "/home/jp/eclipse-workspace/database/MSigDB/c5.go.v7.2.symbols.gmt"

#file.gmt <- "/home/jp/eclipse-workspace/database/MSigDB/c3.tft.gtrd.v7.2.symbols.gmt"
#file.gmt <- "/home/jp/eclipse-workspace/database/MSigDB/c6.all.v7.2.symbols.gmt"
#file.gmt <- "/home/jp/eclipse-workspace/database/MSigDB/jp.gmt"

h.All <- gmtPathways(file.gmt) # 6226
h.All.bis <- read.gmt(file.gmt)

subDir = sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(file.gmt))

i=0
colname <-"pathway"
for (file in files.list){
       
  asbolutepath2file <- glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}/{file}")

  file <- gsub(file,pattern= "-fgsea-gene-ratios.txt",replacement="",perl=T)
  file <- gsub(file,pattern= "JP-",replacement="",perl=T)

  subpath <- file.path(glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}/{file}/"), subDir)
  dir.create(subpath, recursive = TRUE)
  
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
  egmt2 <- GSEA(ranks_decreasing, TERM2GENE=h.All.bis, by = "fgsea",verbose=TRUE ,nPermSimple = 10000 ,minGSSize  = 10, maxGSSize  = 325 , eps = 0,  pvalueCutoff = 1)
  # You dont need the whole object to be written
  write.table(egmt2, file=glue("{subpath}/{file}-full-gsea-clusterprofiler.tsv"),quote=F,row.names=F,sep="\t")

  # Yeah I do it again (I know.Don't say a fucking word moron.)
  fgseaRes     <- fgseaMultilevel(pathways=h.All, stats=ranks,eps=0, nPermSimple = 10000 ,minSize  = 10, maxSize  = 325)
  
  for (pathway in names(h.All)){
    if (pathway %in% c("BLUM_RESPONSE_TO_SALIRASIB_DN","BLUM_RESPONSE_TO_SALIRASIB_UP" )){
      #print(pathway)
      #print(h.All[[pathway]])
      p<- plotEnrichment(h.All[[pathway]], ranks) + labs(title=pathway)
    
      png(file=glue("{subpath}/{file}-{pathway}.png"))
      print(p)
      dev.off()
      
    }  
  }
  
  fgseaResTidy <- fgseaRes %>% as_tibble()

  fwrite( filter(fgseaResTidy ,padj <= 1), file=glue("{subpath}/{file}-full-fgsea.tsv"), sep="\t", sep2=c("", "/", ""))

  topPathwaysUp <- fgseaRes[ES > 0][head(order(pval), n=10), pathway]

  topPathwaysDown <- fgseaRes[ES < 0][head(order(pval), n=10), pathway]

  
  topPathways <- c(topPathwaysUp, rev(topPathwaysDown))
  
  png(file=glue("{subpath}/{file}-top-global.png"),width=900)
  plotGseaTable(h.All[topPathways], ranks, fgseaRes, gseaParam=0.5) 
  dev.off()
  
  
  png(file=glue("{subpath}/{file}-global.png"),width=900)
  print(ggplot( filter(fgseaResTidy ,abs(NES) >= 1.7), aes(reorder(pathway, NES), NES)) +
    geom_col( aes(fill=padj<0.01),color="#5F6368") +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",title=file) + 
    theme_minimal())
  dev.off()
  
  test <- c("PRC2_EZH2_UP.V1_UP","PRC2_EZH2_UP.V1_DN","BLUM_RESPONSE_TO_SALIRASIB_DN","BLUM_RESPONSE_TO_SALIRASIB_UP","E2F1_UP.V1_UP","E2F1_UP.V1_DN")
  
  #write.table(dataframe.expression,glue("/home/jp/eclipse-workspace/Toulouse/data/results/{cellline}/{file}-expresssion.tsv"),quote=F,row.names=T,sep="\t")

  heatplot_up <- heatplot(egmt2,foldChange=ranks, showCategory =   topPathwaysUp)
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

#H4.dormance = H4 DTC vs H4 Untreated
#H4.expansion = H4 DTEC vs H4 DTC Untreated

#PC9.dormance = PC9 DTC vs PC9 Untreated
#PC9.expansion = PC9 DTEC vs PC9 DTC Untreated

#H3.control  = H3 21j vs H3 Untreated

dataset.padj <-  dataset.padj %>% select(-contains(".y"))
dataset.nes  <-  dataset.nes %>% select(-contains(".y"))

dataset.padj <-rename(dataset.padj ,h3.contrast.21d = `H3255_Erlo_21d-H3255_NT`)
dataset.nes <-rename(dataset.nes ,h3.contrast.21d  =  `H3255_Erlo_21d-H3255_NT`)

#dataset.padj <- rename(dataset.padj,h3.contrast.24h = `H3255_Erlo_24h-H3255_NT`)
#dataset.padj <-rename(dataset.padj ,h4.dormance.24h = `H4006_Erlo_24h-H4006_CT_24h`) 
#dataset.nes <-rename(dataset.nes ,h3.contrast.24h  = `H3255_Erlo_24h-H3255_NT`)
#dataset.nes <-rename(dataset.nes ,h4.dormance.24h  =  `H4006_Erlo_24h-H4006_CT_24h`)

if (cellline=="pc9"){
  dataset.padj <-rename(dataset.padj ,pc9.dormance = `PC9_Erlo_DTC-PC9_CT`)
  dataset.padj <-rename(dataset.padj ,pc9.expansion= `PC9_3_DTEC_75d-PC9_Erlo_DTC`)
  dataset.nes <-rename(dataset.nes ,pc9.dormance =  `PC9_Erlo_DTC-PC9_CT`)
  dataset.nes <-rename(dataset.nes ,pc9.expansion =   `PC9_3_DTEC_75d-PC9_Erlo_DTC`)
}

if (cellline=="h4"){
  dataset.padj <-rename(dataset.padj ,h4.dormance = `H4006_Erlo_DTC-H4006_CT_24h`)
  dataset.padj <-rename(dataset.padj ,h4.expansion = `H4006_Erlo_DTEC-H4006_Erlo_DTC`)
  dataset.nes <-rename(dataset.nes ,h4.expansion = `H4006_Erlo_DTEC-H4006_Erlo_DTC`)
  dataset.nes <-rename(dataset.nes ,h4.dormance  =  `H4006_Erlo_DTC-H4006_CT_24h`)
}


write.table(dataset.padj,file=final.file.padj,quote=F,row.names=F,sep="\t")
write.table(dataset.nes,file=final.file.nes,quote=F,row.names=F,sep="\t")
