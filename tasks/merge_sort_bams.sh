#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--runLen )
shift; runLen=$1
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

module load $SAMTOOLSMOD

# performance stats
if [[ $PERFORM = true ]]; then
    echo -e "$(date): merge_sort_bams.sh is mergeing $runLen bams and running on $(hostname)" &>> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.txt &
fi

# merging
echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tmerge_sort_bams.sh-merge\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

eval samtools merge -c -f --threads $THREADS $CWD/$SM/bam/$SM.bam $CWD/$SM/bam/$SM.{1..$runLen}.bam

if [[ $(wc -c <$CWD/$SM/bam/$SM.bam) -ge 1000 ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tmerge_sort_bams.sh-merge\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tmerge_sort_bams.sh-merge\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

# sort
echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tmerge_sort_bams.sh-sort\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
samtools sort --threads $THREADS -m $(( MEM*1000/THREADS ))M -o $CWD/$SM/bam/$SM.sort.bam $CWD/$SM/bam/$SM.bam


if [[ $(wc -c <$CWD/$SM/bam/$SM.bam) -ge 1000 ]]; then
    samtools flagstat -@ $THREADS $CWD/$SM/bam/$SM.sort.bam &>> $CWD/$SM/metrics/$SM.sort.flagstat.txt
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tmerge_sort_bams.sh-sort\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tmerge_sort_bams.sh-sort\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi


