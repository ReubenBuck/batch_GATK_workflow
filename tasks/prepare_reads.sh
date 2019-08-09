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
--workdir )
shift; CWD=$1
;;
--pigz )
shift; PIGZMOD=$1
;;
--threads )
shift; THREADS=$1
;;
--runLen )
shift; runLen=$1
;;
--samtools )
shift; SAMTOOLSMOD=$1
;;
--perform )
shift; PERFORM=$1
;;
--fastqc )
shift; FASTQCMOD=$1
;;
--workdir )
shift; CWD=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

TASK=${SLURM_ARRAY_TASK_ID}

# assign varaibles based on task ids
R1=${R1arr[$TASK]}
R2=${R2arr[$TASK]}
D1=${D1arr[$TASK]}
D2=${D2arr[$TASK]}

echo $CWD

echo -e "$(date)\nVariables in prepare_reads.sh task $TASK have been assigned as,\nR1 is ${R1}\nR2 is ${R2}\nD1 is ${D1}\nD2 is ${D2}\n" > $CWD/$SM/log/$SM.$TASK.prepare_reads.log


echo -e "$(date)\nprepare_reads.$TASK.sh is running on $(hostname)\n" &>> $CWD/$SM/log/$SM.run.log

# here we can start measuring performance stats
if [[ $PERFORM = true ]]; then
    echo -e "prepare_reads.sh is running on $(hostname)" &>>  $CWD/$SM/metrics/perform_prepare_reads_$SM.$TASK.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_prepare_reads_$SM.$TASK.txt &
elif [[ $PERFORM = false ]]; then
    echo -e "$(date)\nPerformance metrics not recorded\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nPerformance var is $PERFORM, requires true/false, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi


module load $SAMTOOLSMOD
module load $FASTQCMOD
module load $PIGZMOD

# need to make log dirs specific for the file
touch $CWD/$SM/log/$SM.$TASK.prepare_reads.log

if [[ $D2 = *".bam" ]]; then
	echo -e "$(date)\nData is storred in unaligned bam format, converting to fastq\n" &>> $CWD/$SM/log/$SM.run.log
	samtools fastq --threads $THREADS -1 $FQDIR/$SM/$R1 -2 $FQDIR/$SM/$R2 $D1/$D2 &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
elif [[ $D2 = *".cram" ]]; then
	echo -e "$(date)\nData is storred in cram format, converting to fastq\n" &>> $CWD/$SM/log/$SM.run.log
	samtools view -@ $THREADS -b -o $CWD/$SM/tmp/${D1/cram/bam} $D1/$D2 &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
	samtools fastq --threads $THREADS -1 $CWD/$SM/fastq/$R1 -2 $FQDIR/$SM/fastq/$R2 $CWD/$SM/tmp/${D1/cram/bam} &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
else
	echo -e "$(date)\nData is likely storred as compressed fastq, uncompressing\n" &>> $CWD/$SM/log/$SM.run.log
	pigz -cd -p $THREADS $D1/$R1.gz > $CWD/$SM/fastq/$R1
	pigz -cd -p $THREADS $D2/$R2.gz > $CWD/$SM/fastq/$R2
fi


# check if file succesfully uncompressed
if [[ -s $CWD/$SM/fastq/$R1 && -s $CWD/$SM/fastq/$R2 ]]; then
	echo -e "$(date)\nUncompressed read pair files found for task no. $TASK, continuing\n" &>> $CWD/$SM/log/$SM.run.log
else
	echo -e "$(date)\nUncompressed read pair files not found for task no. $TASK, exiting\n" &>> $CWD/$SM/log/$SM.run.log
	scancel -n $SM
fi

	# check quality with fastqc
echo -e "$(date)\nFastqc quality checks for task no. $TASK ...\n" &>> $CWD/$SM/log/$SM.run.log
fastqc -o $CWD/$SM/metrics -d $CWD/$SM/tmp -t $THREADS $CWD/$SM/fastq/$R1 $CWD/$SM/fastq/$R2

echo -e "$(date)\nReads are prepared for mapping for task no. $TASK \n" &>> $CWD/$SM/log/$SM.run.log
echo -e "$(date)\nReads are prepared for mapping for task no. $TASK \n" &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
