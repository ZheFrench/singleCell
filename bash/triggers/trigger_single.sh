#!/bin/bash


########################################################
############     Trigger Filtering      ################
########################################################

/data/villemin/code/singleCell/bash/single.sh /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/CTL/outs CTL
/data/villemin/code/singleCell/bash/single.sh /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/OSI/outs OSI
/data/villemin/code/singleCell/bash/single.sh /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/OSI_TIPI/outs OSI_TIPI
/data/villemin/code/singleCell/bash/single.sh /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/TIPI/outs TIPI


/data/villemin/code/singleCell/bash/align10X.sh
