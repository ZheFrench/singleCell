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

#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/PC9_verte /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/PC9_rouge /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_rouge /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_verte /data/villemin/data/toulouse/scRNA-ALL/CellRanger/

#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Vertes /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Rouges /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Rouges /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Vertes /data/villemin/data/toulouse/scRNA-ALL/CellRanger/

##################################################################
############     qualityControl Bash         #####################
##################################################################

#Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -f TRUE 4006_rouge 4006_verte CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes PC9_rouge PC9_verte

Rscript /data/villemin/code/singleCell/R/qualityControl.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -f TRUE Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI

##################################################################
############     qualityControl Bash         #####################
##################################################################

#Rscript /data/villemin/code/singleCell/R/filterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -f TRUE 4006_rouge 4006_verte CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes PC9_rouge PC9_verte -i /data/villemin/data/toulouse/scRNA-ALL/CellRanger/thresholds.csv

#Rscript /data/villemin/code/singleCell/R/filterCells.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -f TRUE Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI OSI_TIPI_Vertes PC9_rouge PC9_verte -i /data/villemin/data/toulouse/scRNA-ALL/CellRanger/thresholds.csv