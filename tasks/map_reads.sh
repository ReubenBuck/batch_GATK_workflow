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