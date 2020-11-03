ref=/lustre/user/liclab/publicData/igenomes/Homo_sapiens/UCSC/hg19/Sequence/Bowtie2Index/genome
#picard=/home/zhangc/tools/picard-tools-1.119/
picard=/lustre/user/liclab/biotools/picard-tools-1.118/

for i in `ls ../data_1022/*_1.clean.fq.gz`
do
{
	b=`basename  $i`
	sample=`echo $b | cut -f 1 -d '_'`
	bowtie2 -p 2 -t -q -N 1 -L 25 -X 2000 --no-mixed --no-discordant -x $ref -1 $i -2 ${i/1.clean.fq.gz/2.clean.fq.gz} -S $sample".sam" > $sample".bowtie2.log"
	
	samtools view -b -h -o $sample".bam" -q 1 -@ 3 $sample".sam"

	samtools sort -m 2G -@ 3 -o $sample".sort.bam" $sample".bam" 

	java -Xms10g -Xmx10g -XX:ParallelGCThreads=3  -jar $picard/MarkDuplicates.jar TMP_DIR=./ INPUT=$sample".sort.bam" 	\
		OUTPUT=$sample".sort.rmdup.bam" METRICS_FILE=$sample".rmdup.metrics"	\
		VALIDATION_STRINGENCY=SILENT ASSUME_SORTED=true REMOVE_DUPLICATES=true MAX_RECORDS_IN_RAM=50000000

}&
done




