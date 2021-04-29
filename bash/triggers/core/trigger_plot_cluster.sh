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
############      PLOT CLUSTER    ################################
##################################################################

#################  CELLS ################# 

# CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte PC9_rouge PC9_verte 
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes PC9_rouge PC9_verte 4006_rouge 4006_verte

# CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes 4006_rouge 4006_verte

# CTL ROUGE vs CTL VERTE
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte CTL_Rouges CTL_Vertes 
# CTL_Rouges OSI_TIPI_Rouges 
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte CTL_Rouges OSI_TIPI_Rouges 
# CTL_Vertes OSI_TIPI_Vertes
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte CTL_Vertes OSI_TIPI_Vertes 
# OSI_TIPI_Rouges vs OSI_TIPI_Vertes
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte OSI_TIPI_Rouges OSI_TIPI_Vertes 
# 4006_rouge 4006_verte
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte 4006_rouge 4006_verte 

# CTL_Rouges OSI_TIPI_Rouges 
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte CTL_Rouges 4006_rouge 
# CTL_Vertes OSI_TIPI_Vertes
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte CTL_Vertes 4006_verte

# CTL_Rouges OSI_TIPI_Rouges 
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte 4006_rouge OSI_TIPI_Rouges 
# CTL_Vertes OSI_TIPI_Vertes
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.4006_rouge.4006_verte 4006_verte OSI_TIPI_Vertes

#################  PDX #################

# Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_CTL_GRCh38 Human_OSI Human_OSI_TIPI Human_TIPI

# Human_CTL_GRCh38 Human_OSI_TIPI
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI_TIPI Human_CTL_GRCh38 

# Human_OSI Human_OSI_TIPI
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI_TIPI Human_OSI 

# Human_TIPI Human_OSI_TIPI
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI_TIPI Human_TIPI

# Human_OSI Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI Human_CTL_GRCh38 

# Human_TIPI Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/plotCluster.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_TIPI Human_CTL_GRCh38

