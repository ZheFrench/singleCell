
#####################################################################
############     Xenome	(10X Genomicss)					############# 
##################################################################### 

DIR=/data/villemin/genome/xenome/10X
cd ${DIR}

echo $PWD

# ---- Create Symbolique link From ref-----
#ln -s /share/apps/STAR/indexes/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa .
#ln -s /share/apps/STAR/indexes/Mus_musculus/Mus_musculus.GRCm38.89.gtf .

#ln -s /share/apps/STAR/indexes/Homo_sapiens/Homo_sapiens.GRCh38.dna.primary_assembly.fa 
#ln -s /share/apps/STAR/indexes/Homo_sapiens/Homo_sapiens.GRCh38.83.gtf .


####  Cell Ranger : I installed that before but i don't think it's mandatory ####
#conda install -c conda-forge rust
#conda install -c conda-forge go
#conda install -c anaconda bzip2
#conda install -c anaconda make
#conda install -c anaconda clang

#curl -o cellranger-6.0.0.tar.gz "https://cf.10xgenomics.com/releases/cell-exp/cellranger-6.0.0.tar.gz?Expires=1617224154&Policy=eyJTdGF0ZW1lbnQiOlt7IlJlc291cmNlIjoiaHR0cHM6Ly9jZi4xMHhnZW5vbWljcy5jb20vcmVsZWFzZXMvY2VsbC1leHAvY2VsbHJhbmdlci02LjAuMC50YXIuZ3oiLCJDb25kaXRpb24iOnsiRGF0ZUxlc3NUaGFuIjp7IkFXUzpFcG9jaFRpbWUiOjE2MTcyMjQxNTR9fX1dfQ__&Signature=HVcgT4p7TFy~skEsbbb3niptgRKZp14DbeQm982c7SCGOzbCtdbjQrZ2aKY4kVnVWCXXi4MzVToBx6T8Q-U3eiKWBu2oRL7YfDYmOPgtUuNPFZvXfil5YnNQxm4ZpR51kyBdvzqKzO2nha2tzBA45glRSnTasx6086EEelmXRK0PQCIIvzA7AAajk3mfVV-M3gMpQDuvTtX5uJ4I3rZcpXOUOnyj6bZ1LyhmwVJPEjxRVR8RXQ4a-7MJJloU9Spz4GZJvVzHQqlawxkgWq2ljdIBiR3RDYpMOkxCYidDo~A50seMXpKcK3T4dHPOng0jErCwYUHFjoyv5xGvjIp46w__&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA"
#https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation

export PATH=/data/villemin/soft/cellranger-6.0.0:$PATH
# Or use the whole path
#/data/villemin/soft/cellranger-6.0.0/bin/cellranger

#https://github.com/alexdobin/STAR/issues/346


#Human and mouse reference dataset, GRCh38 and mm10

#https://support.10xgenomics.com/single-cell-gene-expression/software/release-notes/build#GRCh38mm10_2020A

#################### SETUP ####################


#References - 3.1.0 (July 24, 2019)

#Human and mouse reference dataset, GRCh38 and mm10 (includes human and mouse V(D)J genes)
# Doc was using GRCh38.93 & GRCm38.93

cellranger mkgtf Homo_sapiens.GRCh38.83.gtf Homo_sapiens.GRCh38.83.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene


cellranger mkgtf Mus_musculus.GRCm38.89.gtf Mus_musculus.GRCm38.89.filtered.gtf \
                 --attribute=gene_biotype:protein_coding \
                 --attribute=gene_biotype:lincRNA \
                 --attribute=gene_biotype:antisense \
                 --attribute=gene_biotype:IG_LV_gene \
                 --attribute=gene_biotype:IG_V_gene \
                 --attribute=gene_biotype:IG_V_pseudogene \
                 --attribute=gene_biotype:IG_D_gene \
                 --attribute=gene_biotype:IG_J_gene \
                 --attribute=gene_biotype:IG_J_pseudogene \
                 --attribute=gene_biotype:IG_C_gene \
                 --attribute=gene_biotype:IG_C_pseudogene \
                 --attribute=gene_biotype:TR_V_gene \
                 --attribute=gene_biotype:TR_V_pseudogene \
                 --attribute=gene_biotype:TR_D_gene \
                 --attribute=gene_biotype:TR_J_gene \
                 --attribute=gene_biotype:TR_J_pseudogene \
                 --attribute=gene_biotype:TR_C_gene


#################### MKREF ####################

#Genome_size_b: 5927818181.0.
#Genome_num_chrs: 260.
#        chr_bin_n_bits = min(18, int(math.log(genome_size_b / genome_num_chrs, 2)))


cellranger mkref --genome=GRCh38 \
				 --nthreads 12 \
                 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh38.83.filtered.gtf \
                 --genome=mm10 \
                 --fasta=Mus_musculus.GRCm38.dna.primary_assembly.fa \
                 --genes=Mus_musculus.GRCm38.89.filtered.gtf \
                 --ref-version=3.1.0


cellranger mkref --genome=GRCh38 \
				 --nthreads 12 \
                 --fasta=Homo_sapiens.GRCh38.dna.primary_assembly.fa \
                 --genes=Homo_sapiens.GRCh38.83.filtered.gtf \
                 --ref-version=3.1.0

# un pti malain avait bouger les liens symboliques.
#ln -s  /share/apps/STAR/indexes/Mus_musculus-GRCm38.89/Mus_musculus.GRCm38.dna.primary_assembly.fa /data/villemin/genome/xenome/10X/Mus_musculus.GRCm38.dna.primary_assembly.fa
#ln -s  /share/apps/STAR/indexes/Mus_musculus-GRCm38.89/Mus_musculus.GRCm38.89.gtf /data/villemin/genome/xenome/10X/Mus_musculus.GRCm38.89.gtf
/data/villemin/soft/cellranger-6.0.0/bin/cellranger mkref --genome=mm10 \
         --nthreads 12 \
                 --fasta=Mus_musculus.GRCm38.dna.primary_assembly.fa \
                 --genes=Mus_musculus.GRCm38.89.filtered.gtf \
                 --ref-version=3.1.0

# Behing the hood for the reference creation
# STAR --runMode genomeGenerate \ 
#--genomeDir /data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10/star \ 
 #--runThreadN 1 \ 
 #--genomeFastaFiles /data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10/fasta/genome.fa \ 
  #--sjdbGTFfile /data/USERS/villemin/genome/xenome/10X/GRCh38_and_mm10/genes/genes.gtf \ 
   #--limitGenomeGenerateRAM 17179869184 \ 
    #--genomeSAsparseD 9 \ 
      #--genomeSAindexNbases 14 \ 
       #--genomeChrBinNbits 18

