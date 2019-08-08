#!/bin/bash

# set defults for unused options
BQSR=false
PERFORM=false

SAMTOOLSMOD=samtools/samtools-1.9-test
GATK=/cluster/software/gatk/gatk-3.8/GenomeAnalysisTK.jar
PICARD=/cluster/software/picard-tools/picard-tools-2.1.1/picard.jar
JAVAMOD=java/openjdk/java-1.8.0-openjdk
PIGZMOD=pigz/pigz-2.4
FASTQCMOD=fastqc/fastqc-0.11.7
BWAMOD=bwa/bwa-0.7.17

# need to read in the config file! 
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
  -c | --config )
    shift; CONFIG=$1
    SM=$(cat $CONFIG | sed '2q;d' | awk {'print $1'}) # pulls out first row of config col 1
    LB=$(echo $(cat $CONFIG | cut -f2) | sed 's/ /,/g') # pulls out entire col 2 and converts to comma sep values
    PL=$(echo $(cat $CONFIG | cut -f3) | sed 's/ /,/g')
    FC=$(echo $(cat $CONFIG | cut -f4) | sed 's/ /,/g')
    LN=$(echo $(cat $CONFIG | cut -f5) | sed 's/ /,/g')
    R1=$(echo $(cat $CONFIG | cut -f6) | sed 's/ /,/g')
    R2=$(echo $(cat $CONFIG | cut -f7) | sed 's/ /,/g')
    D1=$(echo $(cat $CONFIG | cut -f8) | sed 's/ /,/g')
    D2=$(echo $(cat $CONFIG | cut -f9) | sed 's/ /,/g')
    REF=$(cat $CONFIG | sed '2q;d' | awk {'print $10'})
    RECAL=$(cat $CONFIG | sed '2q;d' | awk {'print $11'})
    CWD=$(cat $CONFIG | sed '2q;d' | awk {'print $12'}) # root HPC dir for processing
    BAM=$(cat $CONFIG | sed '2q;d' | awk {'print $13'}) # final destination for files
    GVCF=$(cat $CONFIG | sed '2q;d' | awk {'print $14'})
    METRICS=$(cat $CONFIG | sed '2q;d' | awk {'print $15'})
    LOG=$(cat $CONFIG | sed '2q;d' | awk {'print $16'})
    runLen=$(wc -l $CONFIG | cut -d" " -f1)
    ;;
  -r | --recal )
	shift; BQSR=true
	;;
  -P | --perform )
	shift; PERFORM=true
	;;
  -a | --account )
	shift; ACCOUNT=$1
	;;
  -p | --partition )
	shift; PARTITION=$1
	;;
  -e | --email )
	shift; EMAIL=$1
	;;
  -s | --samtools )
	shift; SAMTOOLSMOD=$1
	;;
  -g | --gatk )
	shift; GATK=$1
	;;
  -t | --picard )
	shift; PICARD=$1
	;;
  -c | --pigz )
	shift; JAVAMOD=$1
	;;
  -f | --fastqc )
shift; FASTQCMOD=$1
	;;
  -b | --bwa )
  shift; BWAMOD=$1
  ;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi


# in performance mode we use the entire node for each task
if [[ $PERFORM = true ]]; then
	EXCLUSIVE="--exclusive"
elif [[ $PERFORM = false ]]; then
	EXCLUSIVE=""
fi

# need to check for partitions
# and account and email


echo $LN
echo $R1
echo $REF
echo $PERFORM

 


# prepare dirs
sbatch \
--mem=10g --time=2-00:00 --nodes=1 --ntasks=1 \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL,CANCEL --output=prepare_dirs-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/prepare_dirs.sh --sample $SM \
--platform $PL --flowcell $FC --lane $LN --library $LB \
--read1 $R1 --read2 $R2 --path1 $D1 --path2 $D2 \
--ref $REF --workdir $CWD --recal $RECAL --bam $BAM \
--gvcf $GVCF --metrics $METRICS --log $LOG \
--bwa $BWAMOD --fastqc $FASTQCMOD --pigz $PIGZMOD --java $JAVAMOD \
--samtools $SAMTOOLSMOD --gatk $GATK --picard $PICARD \
--runLen $runLen --perform $PERFORM --bqsr $BQSR \




# Prepare reads for mapping
#sbatch \
#--mem 10g \
#--time 2-00:00 \
#--job.name=$SM \
#--array=1-$runLen \
#-N1 \
#-n20 \
#prepare_reads.sh --platform=$PL --flowcell=$FL \
#--lane=$LN --library=$LB --sample=$SM \
#--read1=$R1 --read2=$R2 --path1=$D1 --path2=$D2 \
#--ref=$REF --workdir=$CWD --recal=$RECAL --threads 20
 

# Map reads

#sbatch --array=1-$runLen map_reads.sh








