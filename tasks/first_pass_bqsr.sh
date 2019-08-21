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
--recal )
shift; RECAL=$1
;;
--perform )
shift; PERFORM=$1
;;
--threads )
shift; THREADS=$1
;;
--workdir )
shift; CWD=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $JAVAMOD

if [[ $PERFORM = true ]]; then
    echo -e "$(date): first_pass_bqsr.sh is running on $(hostname)" &>>  $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt &
fi

echo -e "$(date)\tbegin\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

# make a list of bams to pass to BQSR
LOCI=$(echo $(ls ${REF%/*}/target_loci) |  sed 's/.intervals//g' | sed 's/ /,/g')
eval ls $CWD/$SM/bam/$SM.{$(echo $LOCI)}.realign.bam > $CWD/$SM/tmp/$SM.bams.list


java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
	-nct $THREADS \
	-T BaseRecalibrator \
	-R $REF \
	-I $CWD/$SM/tmp/$SM.bams.list \
	-knownSites $RECAL \
	-o $CWD/$SM/metrics/$SM.recal_data.table

if [[ -s $CWD/$SM/metrics/$SM.recal_data.table ]]; then
    echo -e "$(date)\tend\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\tfail\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi
