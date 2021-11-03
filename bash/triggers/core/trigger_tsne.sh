#!/bin/bash

#################################################################
#
#date: April 01, 2021
#platform: Ubuntu 18.04
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# trigger_qualityControl.sh
# Usage :
#
# ./trigger_qualityControl.sh 
#
# Description :
#
# Will call Rscript to execute quality controls.
#

#Cell ranger directories need to be in the same dir.

#################  CELLS #################

Rscript /data/villemin/code/singleCell/R/tsne.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger CTL_Rouges CTL_Vertes 4006_rouge 4006_verte 

#################  PDX ################# OSI_TIPI_Rouges OSI_TIPI_Vertes 

#Rscript /data/villemin/code/singleCell/R/tsne.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI


#################  SINGLE CELL  #################
Rscript /data/villemin/code/singleCell/R/scfGSEA.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger/DE/  -a H

Rscript /data/villemin/code/singleCell/R/heatmap.fGSea.R -a H -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger/DE/GSEA/  4006_rouge_CTL_Rouges  4006_verte_CTL_Vertes  
Rscript /data/villemin/code/singleCell/R/heatmap.fGSea.R -a C5 -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger/DE/GSEA/  4006_rouge_CTL_Rouges  4006_verte_CTL_Vertes  -f

#4006_verte_CTL_Vertes	

# fGSEA for h4 (hardcoded)
#Rscript /data/villemin/code/singleCell/R/heatmap.fGSEA.R -a C5

# Extract genes core enriched for specific pathways in all comparisons for one cell line
Rscript /data/villemin/code/singleCell/R/extract.genes.from.Pathway.R -a C5 

# Plot matrice for heatmap for this set of genes singleCell
Rscript /data/villemin/code/singleCell//R/heatmap-customGenes.R -p /data/villemin/code/Toulouse-rnaseq/data/results/h4/C5.subset.genes.txt
Rscript /data/villemin/code/singleCell//R/heatmap-customGenes.R -p /data/villemin/code/Toulouse-rnaseq/data/results/h4/test.39.txt

# Plot matrice for heatmap for this set of genes rnaseq
Rscript /data/villemin/code/Toulouse-rnaseq/src/R/heatmap-customGenes.R -p /data/villemin/code/Toulouse-rnaseq/src/bash/subset.C5.genes.txt

# Plot heatmap of GSEA 
#Rscript /data/villemin/code/Toulouse-rnaseq/src/R/heatmap.rnaseq.fGSea.R -a C5 


exit 0 

#### RNA SEQ #####

# fGSEA
#Rscript /data/villemin/code/Toulouse-rnaseq/src/R/fGSEA.R -a H -c h4

# Extract genes core enriched for specific pathways in all comparisons for one cell line
Rscript /data/villemin/code/Toulouse-rnaseq/src/R/extract.genes.from.Pathway.R -a C5 -c h4

# Plot matrice for heatmap for this set of genes
Rscript /data/villemin/code/Toulouse-rnaseq/src/R/heatmap-customGenes.R -p /data/villemin/code/Toulouse-rnaseq/data/results/h4/C5.subset.genes.txt
Rscript /data/villemin/code/Toulouse-rnaseq/src/R/heatmap-customGenes.R -p /data/villemin/code/Toulouse-rnaseq/data/results/h4/test.39.txt

# Plot heatmap of GSEA 
#Rscript /data/villemin/code/Toulouse-rnaseq/src/R/heatmap.rnaseq.fGSea.R -a C5 
