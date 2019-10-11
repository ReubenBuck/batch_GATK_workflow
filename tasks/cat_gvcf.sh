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
--gvcf )
shift; GVCF=$1
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
--memrequest )
shift; MEM=$1
;;
--array-len )
shift; ARRAYLEN=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $JAVAMOD

TASKS=$(echo $(seq -f "%05g" 1 $ARRAYLEN) | sed 's/ /,/g')

# take two different inputs
if [[ BQSR = true ]]; then
	inStatus=recal
elif [[ BQSR = false ]]; then
	inStatus=realign
fi

if [[ $PERFORM = true ]]; then
    echo -e "$(date) cat_gvcf.sh is running on $(hostname)" &>> $CWD/$SM/metrics/perform_cat_gvcf_$SM.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_cat_gvcf_$SM.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_cat_gvcf_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_cat_gvcf_$SM.txt &
fi

echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tcat_gvcf.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

VAR=$(eval echo -e "I=$CWD/$SM/gvcf/$SM.{$(echo $TASKS)}.g.vcf.gz")

java -Djava.io.tmpdir=$CWD/$SM/tmp -Xmx$(( MEM*1000/100*95 ))M -jar $PICARD SortVcf $VAR O=$CWD/$SM/gvcf/$SM.g.vcf.gz

if [[ -s $CWD/$SM/gvcf/$SM.g.vcf.gz ]]; then
    echo -e "$(date)\t${SLURM_JOB_ID}\tend\tcat_gvcf.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tcat_gvcf.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
    scancel -n ${SM}-unmapped
	scancel -n ${SM}-recal-plots
	scancel -n ${SM}-cat-bams
fi
