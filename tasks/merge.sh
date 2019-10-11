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
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $SAMTOOLSMOD

# performance stats
if [[ $PERFORM = true ]]; then
    echo -e "$(date): merge.sh is mergeing $runLen bams and running on $(hostname)" &>> $CWD/$SM/metrics/perform_merge_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_merge.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_merge.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_merge.txt &
fi

# merging
echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tmerge.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

eval samtools merge -c -f --threads $THREADS $CWD/$SM/bam/$SM.sort.bam $CWD/$SM/bam/$SM.{1..$runLen}.sort.bam

if [[ $(wc -c <$CWD/$SM/bam/$SM.sort.bam) -ge 1000 ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tmerge.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tmerge.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi