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
--exome )
shift; EXOME=$1
;;
--array-len )
shift; ARRAYLEN=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

echo -e "Node: $(hostname)\n"

module load $JAVAMOD

if [[ $PERFORM = true ]]; then
    echo -e "$(date): first_pass_bqsr.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_first_pass_bqsr_$SM.txt &
fi

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

#ADD exome filter
if [[ ! -z $EXOME ]]; then
    DOEXOME=$(echo -e "$EXOME -ip 100")
else
    DOEXOME=$(echo -e "$CWD/$SM/tmp/$SM.bqsr.train.bed")
fi

echo $DOEXOME

# make a list of bams to pass to BQSR
TASKS=$(echo $(seq -f "%05g" 1 $ARRAYLEN) | sed 's/ /,/g')
eval ls $CWD/$SM/bam/$SM.{$(echo $TASKS)}.realign.bam > $CWD/$SM/tmp/$SM.bams.list


java -Djava.io.tmpdir=$CWD/$SM/tmp -Xmx$(( MEM*1000/100*95 ))M -jar $GATK \
	-nct $THREADS \
	-T BaseRecalibrator \
	-R $REF \
	-I $CWD/$SM/tmp/$SM.bams.list \
    -L $DOEXOME \
    -XL $CWD/$SM/tmp/$SM.gaps.bed \
	-knownSites $RECAL \
	-o $CWD/$SM/metrics/$SM.recal_data.table

if [[ -s $CWD/$SM/metrics/$SM.recal_data.table ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tfirst_pass_bqsr.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
    scancel -n ${SM}-recal-plots
    scancel -n ${SM}-cat-bams
fi
