#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--sample )
shift; SM=$1
;;
--read1 )
shift; R1=$1
IFS=', ' read -r -a R1arr <<< "$R1"
;;
--read2 )
shift; R2=$1
IFS=', ' read -r -a R2arr <<< "$R2"
;;
--platform )
shift; PL=$1
IFS=', ' read -r -a PLarr <<< "$PL"
;;
--flowcell )
shift; FC=$1
IFS=', ' read -r -a FCarr <<< "$FC"
;;
--lane )
shift; LN=$1
IFS=', ' read -r -a LNarr <<< "$LN"
;;
--library )
shift; LB=$1
IFS=', ' read -r -a LBarr <<< "$LB"
;;
--threads )
shift; THREADS=$1
;;
--samtools )
shift; SAMTOOLSMOD=$1
;;
--bwa )
shift; BWAMOD=$1
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

module load $SAMTOOLSMOD
module load $BWAMOD

TASK=${SLURM_ARRAY_TASK_ID}

R1=${R1arr[$TASK]}
R2=${R2arr[$TASK]}

FC=${FCarr[$TASK]}
LN=${LNarr[$TASK]}
LB=${LBarr[$TASK]}
PL=${PLarr[$TASK]}

# here we can start measuring performance stats
if [[ $PERFORM = true ]]; then
    echo -e "$(date): map_reads.sh task $TASK is running on $(hostname)" &>>  $CWD/$SM/metrics/perform_map_reads_$SM.$TASK.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_map_reads_$SM.$TASK.txt &
elif [[ $PERFORM = false ]]; then
    echo -e "$(date)\nPerformance metrics not recorded\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nPerformance var is $PERFORM, requires true/false, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi


#set up read groups
RG="@RG\tID:${LB}.${FC}.${LN}\tPU:${LB}.${FC}.${LN}\tSM:${SM}\tPL:${PL}\tLB:${LB}"

echo -e "$(date)\nBegin mapping for $SM task $TASK with read group:\n$RG\n" &>> $CWD/$SM/log/$SM.run.log
	
(bwa mem -M -R $RG -t $THREADS $REF $CWD/$SM/fastq/$R1 $CWD/$SM/fastq/$R2 | samtools view -Sb - > $CWD/$SM/bam/$SM.$TASK.bam) 2> $CWD/$SM/log/$SM.$TASK.bwa.log

if [[ $(wc -c <$CWD/$SM/bam/$SM.$TASK.bam) -ge 1000 ]]; then
	samtools flagstat -@ $THREADS $CWD/$SM/bam/$SM.$TASK.bam &>> $CWD/$SM/metrics/$SM.$TASK.flagstat.txt
	echo -e "$(date)\nMapping for $SM task $TASK with read group:\n$RG\nis complete\n" &>> $CWD/$SM/log/$SM.run.log
else
	echo -e "$(date)\nBam file not found or too small for $SM task $TASK, exiting\n" &>> $CWD/$SM/log/$SM.run.log
	scancel -n $SM
fi




