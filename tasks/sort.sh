#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--threads )
shift; THREADS=$1
;;
--samtools )
shift; SAMTOOLSMOD=$1
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

echo -e "Node: $(hostname)\n"

module load $SAMTOOLSMOD

TASK=${SLURM_ARRAY_TASK_ID}


# performance stats
if [[ $PERFORM = true ]]; then
    echo -e "$(date): sort.sh is mergeing $runLen bams and running on $(hostname)" &>> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.$TASK.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.$TASK.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.$TASK.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.$TASK.txt &
fi


# clean sort dir
if [[ $(ls $CWD/$SM/bam/ | grep "$SM.$TASK.sort.bam.tmp" | wc -l) -gt 0 ]]; then 
	rm $CWD/$SM/bam/$SM.$TASK.sort.bam.tmp* 
fi

# sort
echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tsort.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
samtools sort --threads $THREADS -m $(( MEM*1000/THREADS/100*90 ))M -o $CWD/$SM/bam/$SM.$TASK.sort.bam $CWD/$SM/bam/$SM.$TASK.bam


if [[ $(wc -c <$CWD/$SM/bam/$SM.$TASK.sort.bam) -ge 1000 ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tsort.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tsort.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi
