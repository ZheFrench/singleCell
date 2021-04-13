
export PATH=/data/villemin/soft/cellranger-6.0.0:$PATH

###############################################################################
############     Cell Ranger -   scRNAseqCells-1    ##########################
###############################################################################

/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=4006_rouge \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-1 \
   --sample=4006_rouge \
   --localcores=16 \
   --localmem=64   
            
/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=4006_verte \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-1 \
   --sample=4006_verte \
   --localcores=16 \
   --localmem=64   

   /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=pc9_rouge \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-1 \
   --sample=pc9_rouge \
   --localcores=16 \
   --localmem=64   

   /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=pc9_verte \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-1 \
   --sample=pc9_verte \
   --localcores=16 \
   --localmem=64   

exit 0



###############################################################################
############     Cell Ranger - PDX GRCh38_and_mm10  / GRCh38     #############
###############################################################################
#FASTQ contains sequence base with character other than [ACGTN].: file: "/data/USERS/villemin/data/toulouse/scRNAseqPDX/Human_CTL_S1_L001_R1_001.fastq.gz", line: 8

/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=Human_CTL_GRCh38_and_mm10 \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=Human_CTL \
   --localcores=16 \
   --localmem=64
 
 /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=Human_CTL_GRCh38 \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=Human_CTL \
   --localcores=16 \
   --localmem=64              


 /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=Human_OSI \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=Human_OSI \
   --localcores=16 \
   --localmem=64        

    /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=Human_OSI_TIPI \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=Human_OSI_TIPI \
   --localcores=16 \
   --localmem=64        

    /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=Human_TIPI \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=Human_TIPI \
   --localcores=16 \
   --localmem=64        


exit 0

#####################################################################
############     Cell Ranger - PDX 		GRCh38_and_mm10	############
##################################################################### 

/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=OSI_TIPI \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=OSI_TIPI \
   --localcores=16 \
   --localmem=64   


         
/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=OSI \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=OSI \
   --localcores=16 \
   --localmem=64   

   /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=TIPI \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=TIPI \
   --localcores=16 \
   --localmem=64   

   /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=CTL \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX \
   --sample=CTL \
   --localcores=16 \
   --localmem=64   

exit 0
#####################################################################
############     Cell Ranger  - scRNAseqCells-2     #################
##################################################################### 

/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=OSI_TIPI_Rouges \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-2 \
   --sample=OSI_TIPI_Rouges \
   --localcores=16 \
   --localmem=64   
            
/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=OSI_TIPI_Vertes \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-2 \
   --sample=OSI_TIPI_Vertes \
   --localcores=16 \
   --localmem=64   

   /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=CTL_Vertes \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-2 \
   --sample=CTL_Vertes \
   --localcores=16 \
   --localmem=64   

   /data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=CTL_Rouges \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqCells-2 \
   --sample=CTL_Rouges \
   --localcores=16 \
   --localmem=64   
   
exit 0