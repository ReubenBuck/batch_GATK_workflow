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


java -Djava.io.tmpdir=$CWD/$SM/tmp -jar $GATK \
-nct 25 \
-ERC GVCF \
-T HaplotypeCaller \
-R $REF \
-I $CWD/$SM/bam/$SM.${TARGET%\.intervals}.$inStatus.bam \
-o $CWD/$SM/gvcf/$SM.${TARGET%\.intervals}.g.vcf.gz


