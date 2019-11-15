#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--bam )
shift; BAM=$1
;;
--log )
shift; LOG=$1
;;
--metrics )
shift; METRICS=$1
;;
--gvcf )
shift; GVCF=$1
;;
--workdir )
shift; CWD=$1
;;
--bqsr )
shift; BQSR=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

echo -e "Node: $(hostname)\n"

if [[ $BQSR = true ]]; then
    inStatus=recal
elif [[ $BQSR = false ]]; then
    inStatus=realign
fi

# collect and colate efficiency metrics
IDS=$(cat $CWD/$SM/log/$SM.run.log | cut -f2 | sort | uniq)

for i in $IDS; do 
JOBNAME=$(grep $i $CWD/$SM/log/$SM.run.log | cut -f 4 | sed "s/.sh-.*$/.sh/g" | uniq | head -n1)
echo -e "$i\tJob name\t$JOBNAME" &>> $CWD/$SM/metrics/$SM.eff.tsv
seff $i | sed "s/: /\t/g" > $CWD/$SM/tmp/jobspec.txt
LINES=$(wc -l $CWD/$SM/tmp/jobspec.txt | cut -f1 -d" ")
printf "$i%.0s\n" $(seq 1 $LINES) > $CWD/$SM/tmp/jobid.txt
paste $CWD/$SM/tmp/jobid.txt $CWD/$SM/tmp/jobspec.txt &>> $CWD/$SM/metrics/$SM.eff.tsv
done

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tcp_files.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

# update permissions if old files exist
if [[ -e $BAM/$SM.bam ]]; then
        chmod -R u+w $BAM/$SM.bam*
fi

if [[ -e $BAM/$SM.unmap.bam ]]; then
        chmod -R u+w $BAM/$SM.unmap.bam*
fi

if [[ -e $BAM/$SM.halfmap.bam ]]; then
        chmod -R u+w $BAM/$SM.halfmap.bam*
fi

if [[ -e $GVCF/$SM.g.vcf.gz ]]; then
        chmod -R u+w $GVCF/$SM.g.vcf.gz*
fi


# bams 
echo $(md5sum $CWD/$SM/bam/$SM.$inStatus.bam | cut -f1 -d' ') $SM.bam > $BAM/$SM.bam.md5 &
echo $(md5sum $CWD/$SM/bam/$SM.unmap.bam | cut -f1 -d' ') $SM.unmap.bam > $BAM/$SM.unmap.bam.md5 &
echo $(md5sum $CWD/$SM/bam/$SM.halfmap.bam | cut -f1 -d' ') $SM.halfmap.bam > $BAM/$SM.halfmap.bam.md5 &

cp $CWD/$SM/bam/$SM.$inStatus.bam $BAM/$SM.bam &
cp $CWD/$SM/bam/$SM.$inStatus.bam.bai $BAM/$SM.bam.bai &
cp $CWD/$SM/bam/$SM.unmap.bam $BAM/$SM.unmap.bam &
cp $CWD/$SM/bam/$SM.unmap.bam.bai $BAM/$SM.unmap.bam.bai &
cp $CWD/$SM/bam/$SM.halfmap.bam $BAM/$SM.halfmap.bam &
cp $CWD/$SM/bam/$SM.halfmap.bam.bai $BAM/$SM.halfmap.bam.bai &

# gvcf
echo $(md5sum $CWD/$SM/gvcf/$SM.g.vcf.gz | cut -f1 -d' ') $SM.g.vcf.gz > $GVCF/$SM.g.vcf.gz.md5 &

cp $CWD/$SM/gvcf/$SM.g.vcf.gz $GVCF/$SM.g.vcf.gz &
cp $CWD/$SM/gvcf/$SM.g.vcf.gz.tbi $GVCF/$SM.g.vcf.gz.tbi &

wait



# set up full path md5sums for final check
touch $CWD/$SM/log/$SM.final.check.txt

echo $(cut -f1 -d' ' $BAM/$SM.bam.md5) $BAM/$SM.bam > $CWD/$SM/tmp/$SM.final.bam.md5
echo $(cut -f1 -d' ' $BAM/$SM.unmap.bam.md5) $BAM/$SM.unmap.bam > $CWD/$SM/tmp/$SM.final.unmap.md5
echo $(cut -f1 -d' ' $BAM/$SM.halfmap.bam.md5) $BAM/$SM.halfmap.bam > $CWD/$SM/tmp/$SM.final.halfmap.md5
echo $(cut -f1 -d' ' $GVCF/$SM.g.vcf.gz.md5) $GVCF/$SM.g.vcf.gz > $CWD/$SM/tmp/$SM.final.gvcf.md5

echo "" > $CWD/$SM/log/$SM.final.check.txt

md5sum -c $CWD/$SM/tmp/$SM.final.bam.md5 >> $CWD/$SM/log/$SM.final.check.txt &
md5sum -c $CWD/$SM/tmp/$SM.final.unmap.md5 >> $CWD/$SM/log/$SM.final.check.txt &
md5sum -c $CWD/$SM/tmp/$SM.final.halfmap.md5 >> $CWD/$SM/log/$SM.final.check.txt &
md5sum -c $CWD/$SM/tmp/$SM.final.gvcf.md5 >> $CWD/$SM/log/$SM.final.check.txt &


wait


# make sure all checks passed
if [[ $(grep "OK" $CWD/$SM/log/$SM.final.check.txt | wc -l) = 4 ]]; then
	echo -e "$(date)\t${SLURM_JOB_ID}\tend\tcp_files.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
	chmod ugo-w $BAM/$SM.* $GVCF/$SM.*
else
	echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tcp_files.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
	rm $BAM/$SM.* $GVCF/$SM.*
	scancel $SM
	sleep 10s
	exit
fi

if [[ -d $LOG/$SM ]]; then
	chmod -R u+w $LOG/$SM
	rm -r $LOG/$SM
fi

if [[ -d $METRICS/$SM ]]; then
	chmod -R u+w $METRICS/$SM
	rm -r $METRICS/$SM
fi

cp -r $CWD/$SM/log $LOG/$SM &
cp -r $CWD/$SM/metrics $METRICS/$SM &

wait


chmod -R ugo-w $LOG/$SM $METRICS/$SM


