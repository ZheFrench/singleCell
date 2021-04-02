
export PATH=/data/villemin/soft/cellranger-6.0.0:$PATH



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
#####################################################################
############     Cell Ranger (10X)					#################
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

