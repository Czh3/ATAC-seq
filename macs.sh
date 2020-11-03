# $1: bam
# $2: out

for i in `ls ../align/*.sort.rmdup.bam`
do
{
	b=`basename $i`
	macs14 -t $i -f BAM -g hs -n ${b/.sort.rmdup.bam/} -w  --nolambda --nomodel -S
}&
done
