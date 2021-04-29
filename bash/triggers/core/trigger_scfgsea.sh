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
############    GSEA          ###############
##################################################################

##################  CELLS ################# 

Rscript /data/villemin/code/singleCell/R/scfGSEA.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger/DE/

Rscript /data/villemin/code/singleCell/R/scfGSEA.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/DE/

