#!/bin/bash

########################################################
############     Mapping			####################
########################################################
function mapping () {

# Note pour plus tard
#find /pathTo/condition -name *.sorted.bam  | 
#xargs -I{} bash -c 'java -jar /home/jean-philippe.villemin/bin/picard.2-6/picard.2-6.jar MarkDuplicates 
#I=$0 O=${0/.sorted.bam/.sorted.marked.bam} M=${0/.bam/.marked_dup_metrics.txt} VALIDATION_STRINGENCY=LENIENT' {} \;
#ls -1 | parallel fastqc {} is in this dir

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

mkdir -p /data/villemin/data/toulouse/scRNAseqPDX/testing/IDs
mkdir -p /data/villemin/data/toulouse/scRNAseqPDX/testing/Mapping

cat /data/villemin/code/singleCell/chr.txt  |  while read line 
#cat /data/villemin/genome/xenome/chrNameLength.txt | grep "_M" | while read line 
do
   # do something with $line here
   chrom=$(echo $line | cut -f1 -d" ")
   end=$(echo $line | cut -f2 -d" ")
   bedlikeregion=$(echo "$chrom:0-$end")
   #export $bedlikeregion

   echo $bedlikeregion
   echo "Test"
   echo "Extract IDs"

   nb_cores=$(find ${DIR}/ -name "*.bam" | wc -l)
   echo "Cores : $nb_cores"
   
   find ${DIR}/ -name *out.bam | env_parallel -j $nb_cores 'samtools view {1} $bedlikeregion | cut -f1 >  {.}.$bedlikeregion.id_reads.txt' 
   echo "Remove Reads IDs"
   find ${DIR}/ -name *out.bam | env_parallel -j $nb_cores 'picard FilterSamReads I={1} O={.}.filetered.bam READ_LIST_FILE={.}.$bedlikeregion.id_reads.txt  FILTER=excludeReadList'
   echo "Remove Previous Bam"
   find ${DIR}/ -name *out.bam | env_parallel -j $nb_cores 'rm {1}'
   echo "Rename new Bam to old one"
   find ${DIR}/ -name *.filetered.bam | env_parallel -j $nb_cores 'mv {1} {= s/.filetered//g; =}'
   echo "Move IDs files"
   find ${DIR}/ -name *.$bedlikeregion.id_reads.txt | parallel mv {1} /data/villemin/data/toulouse/scRNAseqPDX/testing/IDs

    #https://sites.google.com/site/wiki4metagenomics/tools/samtools/converting-bam-to-fastq
	# sort paired read alignment .bam file (sort by name -n)
	#samtools sort -n SAMPLE.bam -o SAMPLE_sorted.bam

	# save fastq reads in separate R1 and R2 files
	#samtools fastq -@ 8 SAMPLE_sorted.bam -1 SAMPLE_R1.fastq.gz -i1 SAMPLE_I1.fastq.gz  -2 SAMPLE_R2.fastq.gz -i2 SAMPLE_I2.fastq.gz  -0 /dev/null -s /dev/null -n
	#samtools view -H test.bam | grep @HD

done
}

########################################################
############     Tracks				####################
########################################################
#find ${DIR}/ -name *filetered.bao | env_parallel -j $nb_cores 'samtools index {1} '

# Pull the bam ? make the wig
#/share/apps/STAR/bin/Linux_x86_64/STAR --runMode inputAlignmentsFromBAM --runThreadN 6 --inputBAMfile /data/villemin/data/toulouse/output/H3255_Erlo/H3255.Erlo_21j.rep3/star_output/SORTED.H3255.Erlo_21j.rep3.bam --outFileNamePrefix /data/villemin/data/toulouse/output/H3255_Erlo/H3255.Erlo_21j.rep3/star_output/H3255.Erlo_21j.rep3_ --outWigType wiggle  --outWigStrand Stranded --outWigNorm RPM
#=========> WIG GOES BIGWIG

########################################################
############     Main				####################
########################################################
source `which env_parallel.bash`

DIR=/data/villemin/data/toulouse/scRNAseqPDX/testing

FILES_FASTQ=`ls -1 ${DIR}/*.fastq.gz | grep "R1"`

#mapping

FILES_BAM=`find ${DIR}/ -name *.bam`

filtering

echo "Done"

