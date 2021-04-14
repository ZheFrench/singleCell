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
#################################################################


##################################################################
############     qualityControl Bash         #####################
##################################################################

# First Batch Erlotinib 2 cell lines
Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger -a 4006_verte -b 4006_rouge -p TRUE

# Second Batch OsiTipi + Control

Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger -a CTL_Vertes -b CTL_Rouges -p TRUE
Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger -a OSI_TIPI_Vertes -b OSI_TIPI_Rouges -p TRUE
 
# PDX Human_CTL_GRCh38 Human_OSI Human_TIPI  Human_OSI_TIPI

Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI -b Human_OSI_TIPI -p TRUE
Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_TIPI -b Human_OSI_TIPI -p TRUE
Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_CTL_GRCh38 -b Human_OSI_TIPI -p TRUE
#[1] 33694  9349
# 33694  4549
# 33694  4800
 
# 33694  4708
# 33694  4800
#33694  9508
