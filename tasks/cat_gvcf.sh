#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--picard )
shift; PICARD=$1
;;
--loci )
shift; LOCI=$1
IFS=', ' read -r -a LOCIarr <<< "$(echo ,$LOCI)"
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
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

module load $JAVAMOD

TASK=${SLURM_ARRAY_TASK_ID}
TARGET=${LOCIarr[$TASK]}

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

echo -e "$(date)\tbegin\tcat_gvcf.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log

VAR=$(eval echo -e "I=$CWD/$SM/gvcf/$SM.{$(echo $LOCI)}.g.vcf.gz" | sed "s/.intervals//g")

java -jar $PICARD SortVcf $VAR O=$CWD/$SM/gvcf/$SM.g.vcf.gz

if [[ -s $CWD/$SM/gvcf/$SM.g.vcf.gz ]]; then
    echo -e "$(date)\tend\tcat_gvcf.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\tfail\tcat_gvcf.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi
