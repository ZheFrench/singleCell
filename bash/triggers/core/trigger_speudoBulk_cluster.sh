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
############     SPEUDO-BULK DIFFERENTIAL          ###############
##################################################################

###### --- PDX

Rscript /data/villemin/code/singleCell/R/DE-speudoBulk_cluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI -b Human_CTL_GRCh38 -m MKI67
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk_cluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_TIPI -b Human_CTL_GRCh38  -m MKI67
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk_cluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_CTL_GRCh38  -m MKI67

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk_cluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_TIPI -b Human_OSI  -m MKI67
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk_cluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_OSI  -m MKI67
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk_cluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_TIPI  -m MKI67

