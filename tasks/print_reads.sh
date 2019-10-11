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
--recal )
shift; RECAL=1
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
--threads )
shift; THREADS=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

sleep $((RANDOM % 10))

module load $JAVAMOD

TASK=$(seq -f "%05g" ${SLURM_ARRAY_TASK_ID} ${SLURM_ARRAY_TASK_ID})
TARGET=$SM.$TASK.bed

if [[ $PERFORM = true ]]; then
    echo -e "$(date): first_pass_bqsr.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM_$TASK.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM_$TASK.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM_$TASK.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM_$TASK.txt &
fi

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tprint_read.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -Xmx$(( MEM*1000/100*95 ))M -jar $GATK \
-nct $THREADS \
-T PrintReads \
-R $REF \
-L $CWD/$SM/tmp/split_range/$TARGET \
-I $CWD/$SM/bam/$SM.$TASK.realign.bam \
-BQSR $CWD/$SM/metrics/$SM.recal_data.table \
-o $CWD/$SM/bam/$SM.$TASK.recal.bam


if [[ -s $CWD/$SM/bam/$SM.$TASK.recal.bam ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tprint_read.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tprint_read.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
	scancel -n ${SM}-recal-plots
	scancel -n ${SM}-cat-bams
fi



