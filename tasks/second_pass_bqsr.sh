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
--rversion )
shift; RMOD=$1
;;
--ref )
shift; REF=$1
;;
--threads )
shift; THREADS=$1
;;
--recal )
shift; RECAL=$1
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
--exome )
shift; EXOME=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $JAVAMOD
module load $RMOD

if [[ $PERFORM = true ]]; then
    echo -e "$(date): second_pass_bqsr.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_second_pass_bqsr_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_second_pass_bqsr_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_second_pass_bqsr_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_second_pass_bqsr_$SM.txt &
fi

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tsecond_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log


if [[ ! -z $EXOME ]]; then
    DOEXOME=$(echo -e "$EXOME -ip 100")
else
    DOEXOME=$(echo -e "$CWD/$SM/tmp/$SM.bqsr.test.bed")
fi


java -Djava.io.tmpdir=$CWD/$SM/tmp -Xmx$(( MEM*1000/THREADS/100*95 ))M -jar $GATK \
-nct $THREADS \
-T BaseRecalibrator \
-R $REF \
-I $CWD/$SM/tmp/$SM.bams.list \
-L $DOEXOME \
-XL $CWD/$SM/tmp/$SM.gaps.bed \
-knownSites $RECAL \
-o $CWD/$SM/metrics/$SM.post_recal_data.table \
-BQSR $CWD/$SM/metrics/$SM.recal_data.table

if [[ -s $CWD/$SM/metrics/$SM.post_recal_data.table ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tsecond_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tsecond_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tsecond_pass_bqsr.sh-analyze\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
-T AnalyzeCovariates \
-R $REF \
-before $CWD/$SM/metrics/$SM.recal_data.table \
-after $CWD/$SM/metrics/$SM.post_recal_data.table \
-plots $CWD/$SM/metrics/$SM.recalibration_plots.pdf

if [[ -s $CWD/$SM/metrics/$SM.recalibration_plots.pdf ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tsecond_pass_bqsr.sh-analyze\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tsecond_pass_bqsr.sh-analyze\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi
