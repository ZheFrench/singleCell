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
# ./main.sh 
#
# Description :
#
# Will call other bash scripts that will call Rscripts
#

#Cell ranger directories need to be in the same dir.

#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/PC9_verte /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/PC9_rouge /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_rouge /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_verte /data/villemin/data/toulouse/scRNA-ALL/CellRanger/

#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Vertes /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Rouges /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Rouges /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Vertes /data/villemin/data/toulouse/scRNA-ALL/CellRanger/

##################################################################
############     QC PLOT  FILTER       ###########################
##################################################################
./core/trigger_qualityControl.sh

##################################################################
############     SPEUDO-BULK DIFFERENTIAL          ###############
##################################################################
./core/trigger_speudoBulk.sh

########################################################
############      CLUSTER        #######################
#######################################################
./core/trigger_cluster.sh

#############################################################
############      PLOT CLUSTER        #######################
#############################################################

./core/trigger_plot_cluster.sh

#############################################################
############      PLOT CLUSTER        #######################
#############################################################

./core/trigger_scfGSEA.sh
