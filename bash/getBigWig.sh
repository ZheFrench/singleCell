#!/bin/bash

#################################################################
#
#date: April 01, 2021
#platform: Ubuntu 18.04
# author: Villemin Jean-Philippe
# team: Bioinformatique et biologie des systèmes du cancer : J. Colinge 
# Institute : IRCM
#
# getBigWig.sh
# Usage :
#
# ./getBigWig.sh
#
# Description :
#
# Create BigWig files for ucsc browser.
#
#################################################################


#conda activate Pit-3.7.7
#conda install -c bioconda deeptools

#bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_verte/outs/possorted_genome_bam.bam -o /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_verte/outs/ERLO_4006_verte.bw
bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_rouge/outs/possorted_genome_bam.bam -o /data/villemin/data/toulouse/scRNAseqCells-1/CellRanger/4006_rouge/outs/ERLO_4006_rouge.bw

bamCoverage -p 16 -p 16 -b /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Rouges/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Rouges/outs/CTL_Rouges.bw
bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Vertes/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/CTL_Vertes/outs/CTL_Vertes.bw

bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Rouges/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Rouges/outs/OSI_TIPI_Rouges.bw
bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Vertes/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqCells-2/CellRanger/OSI_TIPI_Vertes/outs/OSI_TIPI_Vertes.bw

# PDX Human_CTL_GRCh38 Human_OSI Human_TIPI  Human_OSI_TIPI

bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_CTL_GRCh38/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_CTL_GRCh38/outs/Human_CTL_GRCh38.bw
bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_OSI/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_OSI/outs/Human_OSI.bw
bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_TIPI/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_TIPI/outs/Human_TIPI.bw
bamCoverage -p 16 -b /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_OSI_TIPI/outs/possorted_genome_bam.bam -o  /data/villemin/data/toulouse/scRNAseqPDX/CellRanger/Human_OSI_TIPI/outs/Human_OSI_TIPI.bw