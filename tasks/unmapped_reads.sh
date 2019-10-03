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
    echo -e "$(date): unmapped_reads.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_unmapped_reads_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_unmapped_reads_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_unmapped_reads_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_unmapped_reads_$SM.txt &
fi


echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tunmapped_reads.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

samtools view -@ $THREADS -Sbh -f 12 -F 0x900 $CWD/$SM/bam/$SM.markdup.bam > $CWD/$SM/bam/$SM.unmap.bam

samtools index -@ $THREADS $CWD/$SM/bam/$SM.unmap.bam


samtools view -@ $THREADS -Sbh -f 4 -F 8 -F 0x900 $CWD/$SM/bam/$SM.markdup.bam > $CWD/$SM/tmp/halfmapped.f4F8.bam
samtools view -@ $THREADS -Sbh -f 8 -F 4 -F 0x900 $CWD/$SM/bam/$SM.markdup.bam > $CWD/$SM/tmp/halfmapped.f8F4.bam
samtools cat $CWD/$SM/tmp/halfmapped.{f4F8,f8F4}.bam | samtools sort -@ $THREADS > $CWD/$SM/bam/$SM.halfmap.bam

samtools index -@ $THREADS $CWD/$SM/bam/$SM.halfmap.bam


if [[ $(wc -c <$CWD/$SM/bam/$SM.unmap.bam) -ge 1000 && $(wc -c <$CWD/$SM/bam/$SM.halfmap.bam) -ge 1000 ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tunmapped_reads.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tunmapped_reads.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
	scancel -n ${SM}-recal-plots
	scancel -n ${SM}-cat-bams
fi

