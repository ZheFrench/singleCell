
#export PATH=/data/villemin/soft/cellranger-6.0.0:$PATH


#####################################################################
############     Cell Ranger (10X)					#################
##################################################################### 

/data/villemin/soft/cellranger-6.0.0/bin/cellranger count --id=PDX \
--transcriptome=/data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10 \
   --fastqs=/data/villemin/data/toulouse/scRNAseqPDX/testing \
   --sample=CTL \
   --localcores=16 \
   --localmem=64   
            
