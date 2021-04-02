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

echo "############  FILTERING  ###############"

mkdir -p ${DIR}/IDs
mkdir -p ${DIR}/Mapping

cat /data/villemin/code/singleCell/chr.txt  |  while read line 
do
   # do something with $line here
   chrom=$(echo $line | cut -f1 -d" ")
   end=$(echo $line | cut -f2 -d" ")
   bedlikeregion=$(echo "$chrom:0-$end")
   
   nb_cores=$(find ${DIR}/ -name "*.bam" | wc -l)
   echo "Cores : $nb_cores"
   
   # Note pour plus tard
   #find /pathTo/condition -name *.sorted.bam  | xargs -I{} bash -c 'java -jar /home/jean-philippe.villemin/bin/picard.2-6/picard.2-6.jar MarkDuplicates 
   #I=$0 O=${0/.sorted.bam/.sorted.marked.bam} M=${0/.bam/.marked_dup_metrics.txt} VALIDATION_STRINGENCY=LENIENT' {} \;

   echo "Extract IDs for $bedlikeregion"
   #find ${DIR}/ -name *bam.bam | env_parallel -j $nb_cores 'samtools view {1} $bedlikeregion | cut -f1 >  {.}.$bedlikeregion.id_reads.txt' 


done

mv ${DIR}/*.id_reads.txt ${DIR}/IDs
cat ${DIR}/IDs/*.id_reads.txt > ${DIR}/IDs/all.reads.txt

echo "Remove Reads IDs - Use 32 Cores Available."
bamfile=$(find ${DIR}/ -name *bam.bam) 
echo $bamfile

#REMOVE READS MAPPED ON  MOUSE
#/data/villemin/anaconda3/envs/Pit-3.7.7/bin/java -XX:ParallelGCThreads=32 -jar /data/villemin/anaconda3/envs/Pit-3.7.7/share/picard-2.25.1-0/picard.jar FilterSamReads I=$bamfile O=${bamfile%%.*}.Filtered.bam READ_LIST_FILE=${DIR}/IDs/all.reads.txt FILTER=excludeReadList

#SORTING  #https://sites.google.com/site/wiki4metagenomics/tools/samtools/converting-bam-to-fastq
echo "Sort Bam By Name - Use 32 Cores Available."

#samtools sort -@ 32  -n ${bamfile%%.*}.Filtered.bam -o ${bamfile%%.*}.Filtered.sortedByName.bam

#EXPORT
echo "Create Pure Human FastQ Reads - Use 32 Cores Available."

samtools fastq -@ 32  ${bamfile%%.*}.Filtered.sortedByName.bam -0 ${bamfile%%.*}_R2.fastq.gz -i1 ${bamfile%%.*}_I1.fastq.gz  -i2 ${bamfile%%.*}_I2.fastq.gz  -s ${bamfile%%.*}_S2.fastq.gz -n


echo "Done"

}

########################################################
############     Tracks				####################
########################################################
#find ${DIR}/ -name *filetered.bao | env_parallel -j $nb_cores 'samtools index {1} 'q

# Pull the bam ? make the wig
#/share/apps/STAR/bin/Linux_x86_64/STAR --runMode inputAlignmentsFromBAM --runThreadN 6 --inputBAMfile /data/villemin/data/toulouse/output/H3255_Erlo/H3255.Erlo_21j.rep3/star_output/SORTED.H3255.Erlo_21j.rep3.bam --outFileNamePrefix /data/villemin/data/toulouse/output/H3255_Erlo/H3255.Erlo_21j.rep3/star_output/H3255.Erlo_21j.rep3_ --outWigType wiggle  --outWigStrand Stranded --outWigNorm RPM
#=========> WIG GOES BIGWIG

########################################################
############     Main				############################
########################################################
source `which env_parallel.bash`

DIR=/data/villemin/data/toulouse/scRNAseqPDX/testing

FILES_FASTQ=`ls -1 ${DIR}/*.fastq.gz | grep "R1"`

#mapping
DIR=$1
 
FILES_BAM=`find ${DIR}/ -name *.bam`

filtering

    
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