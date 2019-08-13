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
shift; RECAL=1
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
    echo -e "$(date): second_pass_bqsr.sh is running on $(hostname)" &>>  $CWD/$SM/metrics/perform_second_pass_bqsr_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_second_pass_bqsr_$SM.txt &
fi

echo -e "$(date)\nSecond pass BQSR for sample $SM\n" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
-nct 20 \
-T BaseRecalibrator \
-R $REF \
-I $CWD/$SM/tmp/$SM.bams.list \
-knownSites $RECAL \
-o $CWD/$SM/metrics/$SM.post_recal_data.table \
-BQSR $CWD/$SM/metrics/$SM.recal_data.table \
--log_to_file $CWD/$SM/log/$SM.bqsr_second.log

if [[ -s $CWD/$SM/bam/$SM.realign.sort.bam.bai ]]; then
    echo -e "$(date)\nSecond pass BQSR for $SM is complete\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nSecond pass BQSR table for $SM is not found or is empty, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi


# this may not work! It might require R to work
echo -e "$(date)\nAnalyze covariates for sample $SM\n" &>> $CWD/$SM/log/$SM.run.log

java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
-T AnalyzeCovariates \ 
-R $REF \  
-before $CWD/$SM/metrics/$SM.recal_data.table \
-after $CWD/$SM/metrics/$SM.post_recal_data.table \
-plots $CWD/$SM/metrics/$SM.recalibration_plots.pdf

if [[ -s $CWD/$SM/metrics/$SM.recalibration_plots.pdf ]]; then
    echo -e "$(date)\nCovariate analysis for $SM is complete\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nCovariate analysis for $SM is not found or is empty, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi
