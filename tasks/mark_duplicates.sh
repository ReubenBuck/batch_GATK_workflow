#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--picard )
shift; PICARD=$1
;;
--java )
shift; JAVAMOD=$1
;;
--perform )
shift; PERFORM=$1
;;
--workdir )
shift; CWD=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $JAVAMOD


if [[ $PERFORM = true ]]; then
    echo -e "$(date): mark_duplicates.sh is running on $(hostname)" &>>  $CWD/$SM/metrics/perform_mark_duplicates_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_mark_duplicates_$SM.txt &
fi


echo -e "$(date)\nMarking duplicates for sample $SM\n" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $PICARD MarkDuplicates \
TMP_DIR=$CWD/$SM/tmp \
INPUT=$CWD/$SM/bam/$SM.sort.bam \
OUTPUT=$CWD/$SM/bam/$SM.markdup.bam \
METRICS_FILE=$CWD/SM/metrics/$SM.markdup.metrics \
OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 &>> $CWD/$SM/log/$SM.markdup.log

if [[ $(wc -c <$CWD/$SM/bam/$SM.markdup.bam) -ge 1000 ]]; then
    echo -e "$(date)\nMark duplicates for $SM is complete\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\n$SM bam file after mark duplicates not found or too small, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi