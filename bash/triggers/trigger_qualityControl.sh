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

#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/PC9_verte /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/PC9_rouge /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_rouge /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_verte /data/villemin/data/toulouse/scRNA-ALL/CellRanger/

#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Vertes /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Rouges /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Rouges /data/villemin/data/toulouse/scRNA-ALL/CellRanger/
#ln -s /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Vertes /data/villemin/data/toulouse/scRNA-ALL/CellRanger/

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

##################################################################
############     SPEUDO-BULK DIFFERENTIAL          ###############
##################################################################

# ----- Cells
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a 4006_rouge -b CTL_Rouges
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a 4006_verte -b CTL_Vertes

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a 4006_rouge -b 4006_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a PC9_rouge -b PC9_verte

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a PC9_verte -b 4006_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a PC9_rouge -b 4006_rouge

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a CTL_Rouges -b CTL_Vertes
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b OSI_TIPI_Vertes

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b CTL_Rouges
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Vertes -b CTL_Vertes

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Vertes -b 4006_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b 4006_rouge

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Vertes -b PC9_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b PC9_rouge


###### --- PDX

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI -b Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_TIPI -b Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_CTL_GRCh38

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_TIPI -b Human_OSI
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_OSI
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_TIPI


###################################################################################
############      CLUSTER      With PC9_rouge PC9_verte    #######################
###################################################################################

#################  CELLS #################

# With PC9_rouge PC9_verte
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte CTL_Rouges CTL_Vertes OSI_TIPI_Rouges OSI_TIPI_Vertes 4006_rouge 4006_verte 

# CTL ROUGE vs CTL VERTE
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte CTL_Rouges CTL_Vertes 
# CTL_Rouges OSI_TIPI_Rouges 
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte CTL_Rouges OSI_TIPI_Rouges 
# CTL_Vertes OSI_TIPI_Vertes
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte CTL_Vertes OSI_TIPI_Vertes 
# OSI_TIPI_Rouges vs OSI_TIPI_Vertes
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte OSI_TIPI_Rouges OSI_TIPI_Vertes 
# 4006_rouge 4006_verte
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte 4006_rouge 4006_verte 
# PC9_rouge PC9_verte
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte PC9_rouge PC9_verte 
# CTL_Rouges 4006_rouge 
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte CTL_Rouges 4006_rouge 
# CTL_Vertes 4006_verte 
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -i  CTL_Rouges.CTL_Vertes.OSI_TIPI_Rouges.OSI_TIPI_Vertes.PC9_rouge.PC9_verte.4006_rouge.4006_verte CTL_Vertes 4006_verte 

#################  PDX #################

# Human_CTL_GRCh38 Human_OSI_TIPI
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI_TIPI Human_CTL_GRCh38 
# Human_OSI Human_OSI_TIPI
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI_TIPI Human_OSI 
# Human_TIPI Human_OSI_TIPI
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI_TIPI Human_TIPI
# Human_OSI Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_OSI Human_CTL_GRCh38 
# Human_TIPI Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/clusterCells.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -i  Human_CTL_GRCh38.Human_OSI.Human_OSI_TIPI.Human_TIPI Human_TIPI Human_CTL_GRCh38

########################################################################################################
############      CLUSTER  - PC9 Removed FOR INTEGRATED ANALYSIS   ################################
########################################################################################################

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


##################################################################
############     SPEUDO-BULK DIFFERENTIAL          ###############
##################################################################

# ----- Cells
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a 4006_rouge -b CTL_Rouges
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a 4006_verte -b CTL_Vertes

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a 4006_rouge -b 4006_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a PC9_rouge -b PC9_verte

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a PC9_verte -b 4006_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a PC9_rouge -b 4006_rouge

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a CTL_Rouges -b CTL_Vertes
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b OSI_TIPI_Vertes

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b CTL_Rouges
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Vertes -b CTL_Vertes

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Vertes -b 4006_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b 4006_rouge

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Vertes -b PC9_verte
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNA-ALL/CellRanger -a OSI_TIPI_Rouges -b PC9_rouge


###### --- PDX

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI -b Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_TIPI -b Human_CTL_GRCh38
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_CTL_GRCh38

#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_TIPI -b Human_OSI
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_OSI
#Rscript /data/villemin/code/singleCell/R/DE-speudoBulk.R -d /data/villemin/data/toulouse/scRNAseqPDX/CellRanger -a Human_OSI_TIPI -b Human_TIPI

