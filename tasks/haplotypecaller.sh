#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--gatk )
shift; GATK=$1
;;
--loci )
shift; LOCI=$1
IFS=', ' read -r -a LOCIarr <<< "$(echo ,$LOCI)"
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

TASK=${SLURM_ARRAY_TASK_ID}
TARGET=${LOCIarr[$TASK]}

# take two different inputs
if [[ BQSR = true ]]; then
	inStatus=recal
elif [[ BQSR = false ]]; then
	inStatus=realign
fi

if [[ $PERFORM = true ]]; then
    echo -e "$(date) haplotypecaller.sh is running on $(hostname)" &>>  $CWD/$SM/metrics/perform_haplotypecaller_$SM_${TARGET%\.intervals}.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_haplotypecaller_$SM_${TARGET%\.intervals}.txt &
fi

echo -e "$(date)\nHaplotypecaller for sample $SM ${TARGET%\.intervals} \n" &>> $CWD/$SM/log/$SM.run.log


java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
-nct 25 \
-ERC GVCF \
-T HaplotypeCaller \
-R $REF \
-I $CWD/$SM/bam/$SM.${TARGET%\.intervals}.$inStatus.bam \
-o $CWD/$SM/gvcf/$SM.${TARGET%\.intervals}.g.vcf.gz

if [[ -s $CWD/$SM/gvcf/$SM.${TARGET%\.intervals}.g.vcf.gz ]]; then
    echo -e "$(date)\nGenotyping for $SM ${TARGET%\.intervals} is complete\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nGenotyping for $SM ${TARGET%\.intervals} is not found or is empty, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi
