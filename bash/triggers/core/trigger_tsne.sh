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

#Rscript /data/villemin/code/singleCell/R/tsne.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes 4006_rouge 4006_verte 

#################  PDX #################

#Rscript /data/villemin/code/singleCell/R/tsne.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI
