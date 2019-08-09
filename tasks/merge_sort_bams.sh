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
    echo -e "$(date): merge_sort_bams.sh is mergeing $runLen bams and running on $(hostname)" &>>  $CWD/$SM/metrics/perform_merge_sort_bams_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_merge_sort_bams_$SM.txt &
fi

# merging
echo -e "$(date)\nMerge $runLen bams for sample $SM\n" &>> $CWD/$SM/log/$SM.run.log

eval samtools merge -c -f --threads $THREADS $CWD/$SM/bam/$SM.bam $CWD/$SM/bam/$SM.{1..$runLen}.bam &>> $CWD/$SM/log/$SM.merge.log

if [[ $(wc -c <$CWD/$SM/bam/$SM.bam) -ge 1000 ]]; then
    echo -e "$(date)\nMerge for $SM is complete\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\n$SM bam file after merge not found or too small, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

# sort
echo -e "$(date)\nSort $SM bam\n" &>> $CWD/$SM/log/$SM.run.log
samtools sort --threads $THREADS -o $CWD/$SM/bam/$SM.sort.bam $CWD/$SM/bam/$SM.bam &>> $CWD/$SM/log/$SM.sort.log


if [[ $(wc -c <$CWD/$SM/bam/$SM.bam) -ge 1000 ]]; then
    samtools flagstat -@ $THREADS $CWD/$SM/bam/$SM.sort.bam &>> $CWD/$SM/metrics/$SM.sort.flagstat.txt
    echo -e "$(date)\nSort for $SM is complete\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\n$SM bam file after sort not found or too small, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi


