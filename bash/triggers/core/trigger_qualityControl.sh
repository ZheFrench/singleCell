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

##################################################################
############     QC PLOT        ##################################
##################################################################

#Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -f TRUE CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes PC9_rouge PC9_verte 4006_rouge 4006_verte

#Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -f TRUE Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI

##################################################################
############     FILTER          #################################
##################################################################

#Rscript /data/villemin/code/singleCell/R/filterCells.R  -i /data/villemin/data/toulouse/scRNA-ALL/CellRanger/thresholds.csv -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -f TRUE CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes PC9_rouge PC9_verte 4006_rouge 4006_verte

#Rscript /data/villemin/code/singleCell/R/filterCells.R  -i /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/thresholds.csv -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -f TRUE Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI
