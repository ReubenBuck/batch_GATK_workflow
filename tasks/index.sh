#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
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
--threads )
shift; THREADS=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $SAMTOOLSMOD


if [[ $PERFORM = true ]]; then
    echo -e "$(date): index.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_index_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_index_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_index_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_index_$SM.txt &
fi


echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tindex.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

samtools index -@ $THREADS $CWD/$SM/bam/$SM.markdup.bam

if [[ -s $CWD/$SM/bam/$SM.markdup.bam.bai ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tindex.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tindex.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
	scancel -n ${SM}-recal-plots
	scancel -n ${SM}-cat-bams
fi
