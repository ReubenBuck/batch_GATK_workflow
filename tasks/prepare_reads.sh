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
--path1 )
shift; D1=$1
IFS=', ' read -r -a D1arr <<< "$D1"
;;
--path2 )
shift; D2=$1
IFS=', ' read -r -a D2arr <<< "$D2"
;;
--pigz )
shift; PIGZMOD=$1
;;
--threads )
shift; THREADS=$1
;;
--samtools )
shift; SAMTOOLSMOD=$1
;;
--perform )
shift; PERFORM=$1
;;
--workdir )
shift; CWD=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

echo -e "Node: $(hostname)\n"

TASK=${SLURM_ARRAY_TASK_ID}

# assign varaibles based on task ids
R1=${R1arr[$TASK]}
R2=${R2arr[$TASK]}
D1=${D1arr[$TASK]}
D2=${D2arr[$TASK]}

sleep $((RANDOM % 10))

echo -e "$(date)\nVariables in prepare_reads.sh task $TASK have been assigned as,\nR1 is ${R1}\nR2 is ${R2}\nD1 is ${D1}\nD2 is ${D2}\n"


echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tprepare_reads.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log

# here we can start measuring performance stats
if [[ $PERFORM = true ]]; then
    echo -e "$(date): prepare_reads.sh task $TASK is running on $(hostname)" &>> $CWD/$SM/metrics/perform_prepare_reads_$SM.$TASK.txt
    scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_prepare_reads_$SM.$TASK.txt
    echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_prepare_reads_$SM.$TASK.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_prepare_reads_$SM.$TASK.txt &
fi


module load $SAMTOOLSMOD
module load $PIGZMOD

mkdir $CWD/$SM/fastq/task_$TASK

if [[ $D2 = *".bam" ]]; then
	echo -e "$(date)\nData is storred in unaligned bam format, converting to fastq\n"
	samtools fastq --threads $THREADS -1 $CWD/$SM/fastq/task_$TASK/$R1 -2 $CWD/$SM/fastq/task_$TASK/$R2 $D1/$D2
elif [[ $D2 = *".cram" ]]; then
	mkdir $CWD/$SM/tmp/task_$TASK
	echo -e "$(date)\nData is storred in cram format, converting to fastq\n"
	samtools view -@ $THREADS -b -o $CWD/$SM/tmp/task_$TASK/${D2/cram/bam} $D1/$D2 
	samtools fastq --threads $THREADS -1 $CWD/$SM/fastq/task_$TASK/$R1 -2 $CWD/$SM/fastq/task_$TASK/$R2 $CWD/$SM/tmp/${D2/cram/bam}
else
	echo -e "$(date)\nData is likely storred as compressed fastq, uncompressing\n"
	pigz -cd -p $THREADS $D1/$R1.gz > $CWD/$SM/fastq/task_$TASK/$R1
	pigz -cd -p $THREADS $D2/$R2.gz > $CWD/$SM/fastq/task_$TASK/$R2
fi


# check if file succesfully uncompressed
if [[ -s $CWD/$SM/fastq/task_$TASK/$R1 && -s $CWD/$SM/fastq/task_$TASK/$R2 ]]; then
	echo -e "$(date)\nUncompressed read pair files found for task no. $TASK, continuing\n"
	echo -e "$(date)\t${SLURM_JOB_ID}\tend\tprepare_reads.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
else
	echo -e "$(date)\nUncompressed read pair files not found for task no. $TASK, exiting\n"
	echo -e "$(date)\t${SLURM_JOB_ID}\tfail\tprepare_reads.sh\t$SM\t$TASK" &>> $CWD/$SM/log/$SM.run.log
	scancel -n $SM
	scancel -n ${SM}-unmapped
	scancel -n ${SM}-recal-plots
	scancel -n ${SM}-cat-bams
fi




