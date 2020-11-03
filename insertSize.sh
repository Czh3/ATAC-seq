picard=/lustre/user/liclab/biotools/picard-tools-1.118/

for i in `ls *rmdup.bam`
do
{
	java -Xmx8g -jar $picard/CollectInsertSizeMetrics.jar I=$i 	\
		O=$i".insert_size_metrics.txt" H=$i".insert_size_histogram.pdf" M=0.5
}&
done
