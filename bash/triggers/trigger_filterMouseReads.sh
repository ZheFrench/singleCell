#!/bin/bash

#################################################################
#
#date: April 01, 2021
#platform: Ubuntu 18.04
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# trigger_filterMouseReads.sh
# Usage :
#
# trigger_filterMouseReads.sh directory nameCondition
#
# Description :
#
# Will remove mouse reads from PDX and recreate Fastq with pure human reads in order to be aligned by Cell Ranger again.
#
#################################################################

#########################################################
############     Triggers Filtering      ################
#########################################################

/data/villemin/code/singleCell/bash/filterMouseReads.sh /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/CTL/outs CTL
/data/villemin/code/singleCell/bash/filterMouseReads.sh  /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/OSI/outs OSI
/data/villemin/code/singleCell/bash/filterMouseReads.sh  /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/OSI_TIPI/outs OSI_TIPI
/data/villemin/code/singleCell/bash/filterMouseReads.sh  /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/TIPI/outs TIPI


##########################################################
############     Triggers Alignement      ################
##########################################################

/data/villemin/code/singleCell/bash/align10X.sh
