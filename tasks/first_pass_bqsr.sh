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
--memrequest )
shift; MEM=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $JAVAMOD

if [[ $PERFORM = true ]]; then
    echo -e "$(date): first_pass_bqsr.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt &
fi

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

#ADD exome filter
if [[ ! -z EXOME ]]; then
    DOEXOME=$(echo -e "-L $EXOME --interval-set-rule INTERSECTION --interval_padding 100")
fi


# make a list of bams to pass to BQSR
LOCI=$(echo $(ls ${REF%/*}/target_loci) |  sed 's/.intervals//g' | sed 's/ /,/g')
eval ls $CWD/$SM/bam/$SM.{$(echo $LOCI)}.realign.bam > $CWD/$SM/tmp/$SM.bams.list


java -Djava.io.tmpdir=$CWD/$SM/tmp -Xmx${MEM}G -jar $GATK \
	-nct $THREADS \
	-T BaseRecalibrator \
	-R $REF \
	-I $CWD/$SM/tmp/$SM.bams.list \
    -L $CWD/$SM/tmp/$SM.bqsr.train.bed $DOEXOME \
    -XL $CWD/$SM/tmp/$SM.gaps.bed \
	-knownSites $RECAL \
	-o $CWD/$SM/metrics/$SM.recal_data.table \
    -PF $CWD/$SM/metrics/PF.BaseRecalibrator.$SM.txt

if [[ -s $CWD/$SM/metrics/$SM.recal_data.table ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi
