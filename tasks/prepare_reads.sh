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
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi

TASK=${SLURM_ARRAY_TASK_ID}

echo -e "$(date)\nprepare_reads.$TASK.sh is running on $(hostname)\n" &>> $CWD/$SM/log/$SM.run.log

# here we can start measuring performance stats
echo -e "\n\n$(date)\nChecks for performance\n" &>> $CWD/$SM/log/$SM.run.log
if [[ $PERFORM = true ]]; then
    echo -e "$(date)\nSetting up task performance metrics\n" &>> $CWD/$SM/log/$SM.run.log
    echo -e "prepare_reads.sh is running on $(hostname)" >  $CWD/$SM/log/perform_prepare_reads_$SM.$TASK.txt
    vmstat -twn -S m 1 >> $CWD/$SM/log/perform_prepare_dirs_$SM.$TASK.txt &
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

if [[ $D2 = *".bam" ]]; 
	then
	echo -e "$(date)\nData is storred in unaligned bam format, converting to fastq\n" &>> $CWD/$SM/log/$SM.run.log
	samtools fastq --threads $THREADS -1 $FQDIR/$SM/$R1 -2 $FQDIR/$SM/$R2 $D1/$D2 &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
elif [[ $D2 = *".cram" ]]; then
	echo -e "$(date)\nData is storred in cram format, converting to fastq\n" &>> $CWD/$SM/log/$SM.run.log
	samtools view -@ $THREADS -b -o $CWD/$SM/tmp/${D1/cram/bam} $D1/$D2 &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
	samtools fastq --threads $THREADS -1 $CWD/$SM/fastq/$R1 -2 $FQDIR/$SM/fastq/$R2 $CWD/$SM/tmp/${D1/cram/bam} &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
else
	echo -e "$(date)\nData is likely storred as compressed fastq, uncompressing\n" &>> $CWD/$SM/log/$SM.run.log
	pigz -cd -p $THREADS $D1/$R1.gz > $CWD/$SM/fastq/$R1 &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
	pigz -cd -p $THREADS $D2/$R2.gz > $CWD/$SM/fastq/$R2 &>> $CWD/$SM/log/$SM.$TASK.prepare_reads.log
fi




	# check if file succesfully uncompressed
	if [[ -s $FQDIR/$SM/$R1 && -s $FQDIR/$SM/$R2 ]]; 
		then
			echo uncompressed read pair files found &>> $LOGDIR/$START/$SM/$SM.run.log
		else
			echo one or both uncompressed files are not found or are empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
			exit
		fi

	# check quality with fastqc
	echo running fastqc quality checks ... &>> $LOGDIR/$START/$SM/$SM.run.log
	fastqc -o $QUALDIR/$SM $FQDIR/$SM/$R1 $FQDIR/$SM/$R2

	echo begin mapping &>> $LOGDIR/$START/$SM/$SM.run.log
	# perform mapping
	(bwa mem -M -R $RG -t $THREADS $IDX $FQDIR/$SM/$R1 $FQDIR/$SM/$R2 | samtools view -Sb - > $MAPDIR/$SM/$SM.$ROW.bam) 2> $LOGDIR/$START/$SM/$SM.$ROW.aln.log

	#check for bam files and remove fastq files once reads are mapped
	if [ -s $MAPDIR/$SM/$SM.$ROW.bam ]
	then 
		echo bam file found, removing $FQDIR/$SM/{$R1,$R2} &>> $LOGDIR/$START/$SM/$SM.run.log
		rm -r $FQDIR/$SM/{$R1,$R2}
	else
		echo bam file not found or is empty, exiting &>> $LOGDIR/$START/$SM/$SM.run.log
		exit
	fi

	echo end mapping &>> $LOGDIR/$START/$SM/$SM.run.log

done