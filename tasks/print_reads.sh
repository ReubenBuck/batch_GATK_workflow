#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--gatk )
shift; GATK=$1
;;
--loci )
shift; LOCI=$1
IFS=', ' read -r -a LOCIarr <<< "$(echo ,$LOCI)"
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
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

sleep $((RANDOM % 10))

module load $JAVAMOD

TASK=${SLURM_ARRAY_TASK_ID}
TARGET=${LOCIarr[$TASK]}

if [[ $PERFORM = true ]]; then
    echo -e "$(date): first_pass_bqsr.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt &
fi

echo -e "$(date)\tbegin\tprint_read.sh\t$SM\t${TARGET%\.intervals}" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -Xmx${MEM}G -jar $GATK \
-nct 20 \
-T PrintReads \
-R $REF \
-L ${REF%/*}/target_loci/$TARGET \
-I $CWD/$SM/bam/$SM.${TARGET%\.intervals}.realign.bam \
-BQSR $CWD/$SM/metrics/$SM.recal_data.table \
-o $CWD/$SM/bam/$SM.${TARGET%\.intervals}.recal.bam


if [[ -s $CWD/$SM/bam/$SM.${TARGET%\.intervals}.recal.bam ]]; then
    echo -e "$(date)\tend\tprint_read.sh\t$SM\t${TARGET%\.intervals}" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\tfail\tprint_read.sh\t$SM\t${TARGET%\.intervals}" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi



