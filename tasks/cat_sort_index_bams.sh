#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--samtools )
shift; SAMTOOLSMOD=$1
;;
--picard )
shift; PICARD=$1
;;
--rversion )
shift; RMOD=$1
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
--bqsr )
shift; BQSR=$1
;;
--threads )
shift; THREADS=$1
;;
--memrequest )
shift; MEM=$1
;;
--array-len )
shift; ARRAYLEN=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $SAMTOOLSMOD
module load $RMOD

if [[ $BQSR = true ]]; then
    inStatus=recal
elif [[ $BQSR = false ]]; then
    inStatus=realign
fi


if [[ $PERFORM = true ]]; then
    echo -e "$(date): cat_sort_index_bams.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_cat_sort_index_bams_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_cat_sort_index_bams_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_cat_sort_index_bams_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_cat_sort_index_bams_$SM.txt &
fi

TASKS=$(echo $(seq -f "%05g" 1 $ARRAYLEN) | sed 's/ /,/g')

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tcat_sort_index_bams.sh-concat\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

eval samtools cat -o $CWD/$SM/bam/$SM.$inStatus.cat.bam $CWD/$SM/bam/$SM.{$(echo $TASKS)}.$inStatus.bam

if [[ -s $CWD/$SM/bam/$SM.$inStatus.cat.bam ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tcat_sort_index_bams.sh-concat\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tcat_sort_index_bams.sh-concat\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi


echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tcat_sort_index_bams.sh-sort\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

samtools sort -m $(( MEM*1000/THREADS ))M --threads $THREADS -o $CWD/$SM/bam/$SM.$inStatus.bam $CWD/$SM/bam/$SM.$inStatus.cat.bam

if [[ -s $CWD/$SM/bam/$SM.$inStatus.bam ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tcat_sort_index_bams.sh-sort\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tcat_sort_index_bams.sh-sort\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi


echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tcat_sort_index_bams.sh-index\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

samtools index -@ $THREADS $CWD/$SM/bam/$SM.$inStatus.bam

if [[ -s $CWD/$SM/bam/$SM.$inStatus.bam.bai ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tcat_sort_index_bams.sh-index\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tcat_sort_index_bams.sh-index\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi



echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tcat_sort_index_bams.sh-metrics\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $PICARD CollectMultipleMetrics \
TMP_DIR=$CWD/$SM/tmp \
I=$CWD/$SM/bam/$SM.$inStatus.bam \
O=$CWD/$SM/metrics/$SM.multiple_metrics \
R=$REF


if [[ -s $CWD/$SM/metrics/$SM.multiple_metrics.alignment_summary_metrics ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tcat_sort_index_bams.sh-metrics\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tcat_sort_index_bams.sh-metrics\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi