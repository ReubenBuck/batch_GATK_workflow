#!/bin/bash
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
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
--ref )
shift; REF=$1
;;
--workdir )
shift; CWD=$1
;;
--recal )
shift; RECAL=$1
;;
--bam )
shift; BAM=$1
;;
--log )
shift; LOG=$1
;;
--metrics )
shift; METRICS=$1
;;
--gvcf )
shift; GVCF=$1
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
--gatk )
shift; GATK=$1
;;
--picard )
shift; PICARD=$1
;;
--pigz )
shift; PIGZMOD=$1
;;
--bwa )
shift; BWAMOD=$1
;;
--java )
shift; JAVAMOD=$1
;;
--rversion )
shift; RMOD=$1
;;
--bedtools )
shift; BEDTOOLSMOD=$1
;;
--perform )
shift; PERFORM=$1
;;
--bqsr )
shift; BQSR=$1
;;
--taskdir )
shift; TASKDIR=$1
;;
--array-len )
shift; ARRAYLEN=$1
;;
--gap-size )
shift; GAPSIZE=$1
;;
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi



#clean and establish the working dir
if [[ ! -z $(ls $CWD/$SM/) ]]; then
    if [[ $CWD = "" || $SM = "" ]]; then
        echo "Working dir names are empty, exiting"
        scancel -n $SM
    else
        rm -r $CWD/$SM/*
        mkdir -p $CWD/$SM/log
    fi
    echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tprepare_dirs.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
else
    mkdir -p $CWD/$SM/log
    echo -e "$(date)\t${SLURM_JOB_ID}\tbegin\tprepare_dirs.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
fi
mkdir -p $CWD/$SM/fastq
mkdir -p $CWD/$SM/bam
mkdir -p $CWD/$SM/metrics
mkdir -p $CWD/$SM/gvcf
mkdir -p $CWD/$SM/tmp


echo ${SLURM_JOB_ID}

echo -e "prepare_dirs.sh is running on $(hostname)" &>> $CWD/$SM/log/prepare_dirs-${SM}.out

# here we can start measuring performance stats
echo -e "\n\n$(date)\nChecks for performance\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
if [[ $PERFORM = true ]]; then
	echo -e "$(date)\nSetting up task performance metrics\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
	echo -e "$(date): prepare_dirs.sh is running on $(hostname)\n\n\n" > $CWD/$SM/metrics/perform_prepare_dirs_$SM.txt
	scontrol show jobid -dd ${SLURM_JOB_ID} &>> $CWD/$SM/metrics/perform_prepare_dirs_$SM.txt
	echo -e "\n\n\n" &>> $CWD/$SM/metrics/perform_prepare_dirs_$SM.txt
	vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_prepare_dirs_$SM.txt &
elif [[ $PERFORM = false ]]; then
    echo -e "$(date)\nPerformance metrics not recorded\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nPerformance var is $PERFORM, requires true/false, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

# Check for programs

echo -e "\n\n$(date)\nCheck for programs\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out

module load $JAVAMOD
module load $SAMTOOLSMOD
module load $PIGZMOD
module load $BWAMOD
module load $RMOD
module load $BEDTOOLSMOD

Rscript --version; Rexit=$?
if [[ $Rexit = 0 ]]; then
    echo -e "$(date)\nUsing R version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    Rscript --version &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nR did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

bedtools --version; Bedexit=$?
if [[ $Bedexit = 0 ]]; then
    echo -e "$(date)\nUsing bedtools version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    bedtools --version &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nbedtools did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi


java -version; javExit=$?
if [[ $javExit = 0 ]]; then
    echo -e "$(date)\nUsing java version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    java -version &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nJava did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

samtools --version; samExit=$?
if [[ $samExit = 0 ]]; then
    echo -e "$(date)\nUsing samtools version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    samtools --version &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nSamtools did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi


pigz -V; pigExit=$?
if [[ $pigExit = 0 ]]; then
    echo -e "$(date)\nUsing pigz version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    pigz -V &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\Pigz did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

java -jar $GATK -version; gatExit=$?
if [[ $gatExit = 0 ]]; then
    echo -e "$(date)\nUsing GATK version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    java -jar $GATK -version &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\GATK did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

bwa &> $CWD/$SM/tmp/bwa.open.txt
grep Version $CWD/$SM/tmp/bwa.open.txt > $CWD/$SM/tmp/bwa.version.txt
if [[ $(wc -l $CWD/$SM/tmp/bwa.version.txt | cut -d" " -f1) = 1 ]]; then
    echo -e "$(date)\nUsing BWA version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    cat $CWD/$SM/tmp/bwa.version.txt &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nBWA did not report version, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi


java -jar $PICARD MarkDuplicates --version &> $CWD/$SM/tmp/picard.version.txt
if [[ $(wc -l $CWD/$SM/tmp/picard.version.txt | cut -d" " -f1) = 1 ]]; then
    echo -e "$(date)\nUsing picard MarkDuplicates version:" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    cat $CWD/$SM/tmp/picard.version.txt &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nPicard MarkDuplicates did not report version, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

echo -e "\n\n$(date)\nProgram checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out


# check files exist!
# check if unaligned files exist
echo -e "\n\n$(date)\nChecks for unaligned WGS data\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
for i in $(seq 1 $runLen); do
    if [[ ${D2arr[$i]} = *".bam" ]]; then
        if [[ -e ${D1arr[$i]}/${D2arr[$i]} ]]; then
            echo -e "$(date)\nbam file ${D1arr[$i]}/${D2arr[$i]} found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
        else
            echo -e "$(date)\nbam file ${D1arr[$i]}/${D2arr[$i]} not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
            scancel -n $SM
        fi
    elif [[ ${D2arr[$i]} = *".cram" ]]; then
        if [[ -e ${D1arr[$i]}/${D2arr[$i]} ]]; then
            echo -e "$(date)\ncram file ${D1arr[$i]}/${D2arr[$i]} found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
        else
            echo -e "$(date)\ncram file ${D1arr[$i]}/${D2arr[$i]} not found\n", exiting &>> $CWD/$SM/log/prepare_dirs-${SM}.out
            scancel -n $SM
        fi
    else
        if [[ -e ${D1arr[$i]}/${R1arr[$i]}.gz && -e ${D2arr[$i]}/${R2arr[$i]}.gz ]]; then
            echo -e "$(date)\nfastq files ${D1arr[$i]}/${R1arr[$i]}.gz and ${D2arr[$i]}/${R2arr[$i]}.gz found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
        else
            echo -e "$(date)\nfastq files ${D1arr[$i]}/${R1arr[$i]}.gz and ${D2arr[$i]}/${R2arr[$i]}.gz not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
            scancel -n $SM
        fi
    fi
done

echo -e "\n\n$(date)\nData file checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out


# check if ref exists
echo -e "\n\n$(date)\nChecks for reference files\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
if [[ -e $REF ]]; then
    echo -e "$(date)\nref file $REF found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nref file $REF not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

if [[ -e $REF.amb && -e $REF.ann && -e $REF.bwt && -e $REF.pac && -e $REF.sa ]]; then
    echo -e "$(date)\nbwa index files found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nbwa index files not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

if [[ -e ${REF/.fa/}.dict ]]; then
    echo -e "$(date)\ndict file found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\ndict file not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

if [[ -e $REF.fai ]]; then
    echo -e "$(date)\nfai file found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nfai file not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

# check if targets exist
if [[ -d ${REF%/*}/target_loci/ ]]; then
    echo -e "$(date)\ntarget_loci dir found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\ntarget_loci dir not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi 

# check if recal db exists if recal is required
echo -e "\n\n$(date)\nChecks for BQSR\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
if [[ $BQSR = true ]]; then
    if [[ -e $RECAL ]]; then
        echo -e "$(date)\nBQSR DB found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    else
        echo -e "$(date)\nBQSR DB required and not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
        scancel -n $SM
    fi
elif [[ $BQSR = false ]]; then
    echo -e "$(date)\nBQSR is not required from user\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nBQSR var is $BQSR, requires true/false, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

echo -e "\n\n$(date)\nReference file checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out



# check for final destinations and create some tests
# make the dirs for final destinations 
echo -e "\n\n$(date)\nChecks for final destinations\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
if [[ -d $BAM ]]; then
    echo -e "$(date)\nFinal bam dir, $BAM, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nFinal bam dir, $BAM, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

if [[ -d $GVCF ]]; then
    echo -e "$(date)\nFinal GVCF dir, $GVCF, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nFinal GVCF dir, $GVCF, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

if [[ -d $LOG ]]; then
    echo -e "$(date)\nFinal log dir, $LOG, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nFinal log dir, $LOG, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

if [[ -d $METRICS ]]; then
    echo -e "$(date)\nFinal metrics dir, $METRICS, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
else
    echo -e "$(date)\nFinal metrics dir, $METRICS, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out
    scancel -n $SM
fi

echo -e "\n\n$(date)\nFinal destination checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out

echo -e "\n\n$(date)\nPrepare reference intervals" &>> $CWD/$SM/log/prepare_dirs-${SM}.out

# extract gaps from reference index files
Rscript $TASKDIR/process_ref.R $REF > $CWD/$SM/tmp/$SM.gaps.bed
# get chr sizes for bedtools genome file
cut -f1,2 $REF.fai > $CWD/$SM/tmp/$SM.bedtools.genome
# get order of chr names from fai file
cut -f1 $REF.fai > $CWD/$SM/tmp/$SM.names.txt
# create bed file of entire chr lengths
awk '{print $1"\t"0"\t"$2}' $CWD/$SM/tmp/$SM.bedtools.genome > $CWD/$SM/tmp/$SM.genome.bed
# pull out large gap intervals for splitting genome
paste $CWD/$SM/tmp/$SM.gaps.bed <(cat $CWD/$SM/tmp/$SM.gaps.bed | 
    awk -F'\t' '{print $3-$2 }') | awk -v gapsize=$GAPSIZE '$4>=gapsize' | cut -f1-3 > $CWD/$SM/tmp/$SM.large.gaps.bed

# supbtract from genome bed to get list of regions bed
bedtools subtract -a $CWD/$SM/tmp/$SM.genome.bed -b $CWD/$SM/tmp/$SM.large.gaps.bed > $CWD/$SM/tmp/$SM.split.bed

if [[ -s $CWD/$SM/tmp/$SM.gaps.bed && -s $CWD/$SM/tmp/$SM.genome.bed && -s $CWD/$SM/tmp/$SM.large.gaps.bed && -s $CWD/$SM/tmp/$SM.split.bed ]]; then
    echo -e "$(date)\nRef and gap bed files exist"
else
    echo -e "$(date)\nRef and gap bed files not found, exiting"
    scancel $SM
fi


# use region list to split loci into similar sized chunks and sort
mkdir -p $CWD/$SM/tmp/split_range
bedtools split -i $CWD/$SM/tmp/$SM.split.bed -n $ARRAYLEN -p $CWD/$SM/tmp/split_range/$SM
for i in $(ls $CWD/$SM/tmp/split_range/$SM*.bed); do
    bedtools sort -faidx $CWD/$SM/tmp/$SM.names.txt -i $i > $CWD/$SM/tmp/split_range/tmp.bed
    mv $CWD/$SM/tmp/split_range/tmp.bed $i
done
rm $CWD/$SM/tmp/split_range/tmp.bed

if [[ -s $CWD/$SM/tmp/split_range/$SM.00001.bed ]]; then
    echo -e "$(date)\nSplit range found"
else
    echo -e "$(date)\nFirst split range not found, exiting"
    scancel $SM
fi

# get bqsr regions for training and testing
# should speed up analysis by getting a good representation of error profiles
# bqsr train
bedtools random -seed 1 -n 100 -l 1000000 -g $CWD/$SM/tmp/$SM.bedtools.genome | cut -f1-3 | bedtools sort -faidx $CWD/$SM/tmp/$SM.names.txt > $CWD/$SM/tmp/$SM.bqsr.train.bed
# bqsr test
bedtools random -seed 2 -n 100 -l 1000000 -g $CWD/$SM/tmp/$SM.bedtools.genome | cut -f1-3 | bedtools sort -faidx $CWD/$SM/tmp/$SM.names.txt > $CWD/$SM/tmp/$SM.bqsr.test.bed

if [[ -s $CWD/$SM/tmp/$SM.bqsr.train.bed && -s $CWD/$SM/tmp/$SM.bqsr.test.bed ]]; then
    echo -e "$(date)\nTrain and test intervals for bqsr exist"
else
    echo -e "$(date)\nTrain and test intervals not found, exiting"
    scancel $SM
fi



echo -e "\n\n$(date)\nFile and program checks complete, moving on...\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}.out

echo -e "$(date)\t${SLURM_JOB_ID}\tend\tprepare_dirs.sh\t$SM\t" &>> $CWD/$SM/log/$SM.run.log
