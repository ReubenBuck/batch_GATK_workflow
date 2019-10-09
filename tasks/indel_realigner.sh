#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--gatk )
shift; GATK=$1
;;
--java )
shift; JAVAMOD=$1
;;
--ref )
shift; REF=$1
;;
--perform )
shift; PERFORM=$1
;;
--workdir )
shift; CWD=$1
;;
--memrequest )
shift; MEM=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $JAVAMOD

sleep $((RANDOM % 10))

TASK=$(seq -f "%05g" ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_ID})
TARGET=$SM.$TASK.bed

if [[ $PERFORM = true ]]; then
    echo -e "$(date): indel_realigner.sh for $TASK is running on $(hostname)" &>> $CWD/$SM/metrics/perform_indel_realigner_$SM_$TASK.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_indel_realigner_$SM_$TASK.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_indel_realigner_$SM_$TASK.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_indel_realigner_$SM_$TASK.txt &
fi

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tindel_realigner.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
-T IndelRealigner \
-R $REF \
-I $CWD/$SM/bam/$SM.markdup.bam \
-targetIntervals $CWD/$SM/fastq/$SM.indelTarget.intervals \
-L $CWD/$SM/tmp/split_range/$TARGET \
-o $CWD/$SM/bam/$SM.$TASK.realign.bam

if [[ $(wc -c <$CWD/$SM/bam/$SM.$TASK.realign.bam) -ge 1000 ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tindel_realigner.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tindel_realigner.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
	scancel -n ${SM}-recal-plots
	scancel -n ${SM}-cat-bams
fi
