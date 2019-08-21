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
--ref )
shift; REF=$1
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

module load $JAVAMOD
if [[ $PERFORM = true ]]; then
    echo -e "$(date): indel_target_creator.sh is running on $(hostname)" &>>  $CWD/$SM/metrics/perform_indel_target_creator_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_indel_target_creator_$SM.txt &
fi

echo -e "$(date)\tbegin\trealigner_target_creator.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log


java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
-nt $THREADS \
-T RealignerTargetCreator \
-R $REF \
-I $CWD/$SM/bam/$SM.markdup.bam \
-o $CWD/$SM/fastq/$SM.indelTarget.intervals


if [[ $(wc -c <$CWD/$SM/fastq/$SM.indelTarget.intervals) -ge 1000 ]]; then
    echo -e "$(date)\tend\trealigner_target_creator.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\tfail\trealigner_target_creator.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi
