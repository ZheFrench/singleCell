#!/bin/bash

########################################################
############     Mapping			####################
########################################################
function mapping () {

echo "############  MAPPING  ###############"
for fileR1 in ${FILES_FASTQ}
do
	echo ">>> Star Map >>>"
	fileR2="${fileR1//_R1_/_R2_}"
	echo ${fileR1}
	echo ${fileR2}

    name=${fileR1##*/}
    base=${name%.fastq.gz}
    random_dir=$(tr -dc A-Z </dev/urandom | head -c 10 ; echo '')
    echo "/data/villemin/tmp_star/${random_dir}"
    baseUnik="${base//_00*/}"
	echo "${DIR}/${baseUnik}_"

	# Omit 2-pass mapping
	/share/apps/STAR/bin/Linux_x86_64/STAR --runThreadN 6 --genomeDir /data/villemin/genome/xenome/ --readFilesIn ${fileR1} ${fileR2} --outTmpDir "/data/villemin/tmp_star/${random_dir}" --readFilesCommand zcat  --outFileNamePrefix "${DIR}/${baseUnik}_" --outSAMtype BAM SortedByCoordinate --quantMode GeneCounts 

	echo ">>> Samtools Index >>>"

  samtools index "${DIR}/${baseUnik}_Aligned.sortedByCoord.out.bam"

done

}

########################################################
############     Filtering			####################
########################################################
function filtering () {

echo "############  FILTERING Mouse  ###############"

mkdir -p ${DIR}/IDs_MOUSE

cat /data/villemin/code/singleCell/chr_mouse.txt  |  while read line 
do
   # do something with $line here
   chrom=$(echo $line | cut -f1 -d" ")
   end=$(echo $line | cut -f2 -d" ")
   bedlikeregion=$(echo "$chrom:0-$end")
   find ${DIR}/ -name "*.bam" 
   nb_cores=$(find ${DIR}/ -name "*.bam" | wc -l)
   echo "Cores : $nb_cores"
   
   echo "Extract IDs for $bedlikeregion"
   find ${DIR}/ -name *bam.bam | env_parallel -j $nb_cores 'samtools view {1} $bedlikeregion | cut -f1 >  {.}.$bedlikeregion.id_reads.txt' 


done

mv ${DIR}/*.id_reads.txt ${DIR}/IDs_MOUSE
cat ${DIR}/IDs_MOUSE/*.id_reads.txt > ${DIR}/IDs_MOUSE/all.reads.txt


# Finally there is no need of that because I will appy the filtering directly on the initial fastq.

#echo "Remove Reads IDs - Use 32 Cores Available."
#bamfile=$(find ${DIR}/ -name *bam.bam) 
#echo $bamfile


#REMOVE READS MAPPED ON  HUMAN
#/data/villemin/anaconda3/envs/Pit-3.7.7/bin/java -XX:ParallelGCThreads=32 -jar /data/villemin/anaconda3/envs/Pit-3.7.7/share/picard-2.25.1-0/picard.jar FilterSamReads I=$bamfile O=${bamfile%%.*}.Filtered.bam READ_LIST_FILE=${DIR}/IDs_HUMAN/all.reads.txt FILTER=excludeReadList


#SORTING  #https://sites.google.com/site/wiki4metagenomics/tools/samtools/converting-bam-to-fastq
#echo "Sort Bam By Name - Use 32 Cores Available."

#samtools sort -@ 32  -n ${bamfile%%.*}.Filtered.bam -o ${bamfile%%.*}.Filtered.sortedByName.bam

#EXPORT
#echo "Create Pure Human FastQ Reads - Use 32 Cores Available."

#samtools fastq -@ 32  ${bamfile%%.*}.Filtered.sortedByName.bam -0 ${bamfile%%.*}_R2.fastq.gz -n -i 


echo "Done Mouse"

echo "############  FILTERING Human  ###############"

mkdir -p ${DIR}/IDs_HUMAN

cat /data/villemin/code/singleCell/chr_human.txt  |  while read line 
do
   # do something with $line here
   chrom=$(echo $line | cut -f1 -d" ")
   end=$(echo $line | cut -f2 -d" ")
   bedlikeregion=$(echo "$chrom:0-$end")
   
   nb_cores=$(find ${DIR}/ -name "*.bam" | wc -l)
   echo "Cores : $nb_cores"
   
   echo "Extract IDs for $bedlikeregion"
   find ${DIR}/ -name *bam.bam | env_parallel -j $nb_cores 'samtools view {1} $bedlikeregion | cut -f1 >  {.}.$bedlikeregion.id_reads.txt' 


done

mv ${DIR}/*.id_reads.txt ${DIR}/IDs_HUMAN
cat ${DIR}/IDs_HUMAN/*.id_reads.txt > ${DIR}/IDs_HUMAN/all.reads.txt


echo "Done Human"


}


########################################################
############     Main				############################
########################################################
source `which env_parallel.bash`

#DIR=/data/villemin/data/toulouse/scRNAseqPDX/testing

#FILES_FASTQ=`ls -1 ${DIR}/*.fastq.gz | grep "R1"`

# Depreciated
#mapping

DIR=$1
KEYWORD_SAMPLE=$2

BASE_DIR=/data/villemin/data/toulouse/scRNAseqPDX/

echo "Fastq Found for $2"
find /data/villemin/data/toulouse/scRNAseqPDX/ -maxdepth 1 -name "${KEYWORD_SAMPLE}_S*"

filtering

nb_cores=$(find ${BASE_DIR} -maxdepth 1 -name "${KEYWORD_SAMPLE}_S*"  | wc -l)
echo "Cores : $nb_cores"

echo "Mouse FastQ..."
find $BASE_DIR -maxdepth 1 -name "${KEYWORD_SAMPLE}_S*" | env_parallel -j $nb_cores 'seqtk subseq  <(zcat {1}) ${BASE_DIR}/CellRanger/$KEYWORD_SAMPLE/outs/IDs_MOUSE/all.reads.txt | gzip >  $BASE_DIR/Mouse_$(basename {1})'

 echo "Get reads ID unique to Human - This step is not usefull..."
#Find lines only in first file given
comm -23 <(sort ${BASE_DIR}/CellRanger/$KEYWORD_SAMPLE/outs/IDs_HUMAN/all.reads.txt)  <(sort ${BASE_DIR}/CellRanger/$KEYWORD_SAMPLE/outs/IDs_MOUSE/all.reads.txt) > ${BASE_DIR}/CellRanger/$KEYWORD_SAMPLE/IDs.txt

echo "Human FastQ..."
find $BASE_DIR -maxdepth 1 -name "${KEYWORD_SAMPLE}_S*" | env_parallel -j $nb_cores 'seqtk subseq  <(zcat {1}) ${BASE_DIR}/CellRanger/$KEYWORD_SAMPLE/IDs.txt | gzip >  $BASE_DIR/Human_$(basename {1})'

########################################################
############     Tracks       ####################
########################################################
#find ${DIR}/ -name *filetered.bao | env_parallel -j $nb_cores 'samtools index {1} '

# Pull the bam ? make the wig
#/share/apps/STAR/bin/Linux_x86_64/STAR --runMode inputAlignmentsFromBAM --runThreadN 6 --inputBAMfile /data/villemin/data/toulouse/output/H3255_Erlo/H3255.Erlo_21j.rep3/star_output/SORTED.H3255.Erlo_21j.rep3.bam --outFileNamePrefix /data/villemin/data/toulouse/output/H3255_Erlo/H3255.Erlo_21j.rep3/star_output/H3255.Erlo_21j.rep3_ --outWigType wiggle  --outWigStrand Stranded --outWigNorm RPM
#=========> WIG GOES BIGWIG


####################################################################################################################################
############     FASTQ FORMAT       ################################################################################################
####################################################################################################################################

##>>>> R1
#@NB551452:220:HMWJ2BGXH:1:11101:14737:1053 1:N:0:ACCAGACAAC+CCTAGTTCCT
#CGCCAGATCAAACCCACGGCAGGANTAA
#+
#AAAAAEEEEEEEEEEEEEEEEEEE#EEE

##>>>> R2
#@NB551452:220:HMWJ2BGXH:1:11101:14737:1053 2:N:0:ACCAGACAAC+CCTAGTTCCT
#AAATCCACCCCTTANGAGTGCGGCTNCNNCNCNATATNCCCCGCCCGCGTGCCGTTNTNCATAAAATTGTTCTTTGTATGTATTAACTTG
#+
#AAAAAEEEEEEEEE#EEEEEEEAEE#E##E#E#EEAE#EEEEEEEEEEE//EE//E#/#E///<E//E/AE/EE//E/////AA//</A/

#>>>> I1
#@NB551452:220:HMWJ2BGXH:1:11101:14737:1053 1:N:0:ACCAGACAAC+CCTAGTTCCT
#ACCAGACAAC
#+
#AAAAAEEEEE

#>>>> I2
#@NB551452:220:HMWJ2BGXH:1:11101:14737:1053 2:N:0:ACCAGACAAC+CCTAGTTCCT
#CCTAGTTCCT
#+
#A/A6AE/AAA

####################################################################################################################################
############     FASTQ FORMAT       ################################################################################################
####################################################################################################################################


