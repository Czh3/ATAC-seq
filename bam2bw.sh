for i in `ls *sort.rmdup.bam`
do
	bamCoverage -b $i -o ${i/sort.rmdup.bam/bw} -p 5 --normalizeUsing RPKM &
done
