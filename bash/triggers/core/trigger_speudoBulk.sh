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

