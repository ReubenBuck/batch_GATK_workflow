#!/bin/bash
# read in all of the outside flags
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
--platform )
shift; PL=$1
IFS=', ' read -r -a PLarr <<< "$PL"
;;
--flowcell )
shift; FL=$1
IFS=', ' read -r -a FLarr <<< "$FL"
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
esac; shift; done
if [[ "$1" == '--' ]]; then shift; fi


# clean and establish the working dir
echo -e "\n\n$(date)\nChecks for clean working dir\n" &>> $CWD/$SM/log/$SM.run.log
if [[ ! -z $(ls $CWD/$SM/) ]]; then
    rm -r $CWD/$SM/*
    mkdir -p $CWD/$SM/log
    echo -e "$(date)\ndir $CWD/$SM/ is not empty, removing files\n" &>> $CWD/$SM/log/$SM.run.log
else
    mkdir -p $CWD/$SM/log
    echo -e "$(date)\n$dir CWD/$SM/ not found, creating\n" &>> $CWD/$SM/log/$SM.run.log
fi
mkdir -p $CWD/$SM/fastq
mkdir -p $CWD/$SM/bam
mkdir -p $CWD/$SM/metrics
mkdir -p $CWD/$SM/gvcf
mkdir -p $CWD/$SM/tmp

# here we can start measuring performance stats
echo -e "\n\n$(date)\nChecks for performance\n" &>> $CWD/$SM/log/$SM.run.log
if [[ $PERFORM = true ]]; then
    echo -e "$(date)\nSetting up task performance metrics\n" &>> $CWD/$SM/log/$SM.run.log
    vmstat -twn -S m 1 >> $CWD/$SM/log/perform_prepare_dirs_$SM.txt &
elif [[ $PERFORM = false ]]; then
    echo -e "$(date)\nPerformance metrics not recorded\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nPerformance var is $REACL, requires true/false, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

# check files exist!
# check if unaligned files exist
echo -e "\n\n$(date)\nChecks for unaligned WGS data\n" &>> $CWD/$SM/log/$SM.run.log
for i in $(seq 0 $(($runLen - 1))); do
    if [[ $D2 = *".bam" ]]; then
        if [[ -e $D1/$D2 ]]; then
            echo -e "$(date)\nbam file $D2/$D1 found\n" &>> $CWD/$SM/log/$SM.run.log
        else
            echo -e "$(date)\nbam file $D2/$D1 not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
            scancel -n $SM
        fi
    elif [[ $D2 = *".cram" ]]; then
        if [[ -e $D1/$D2 ]]; then
            echo -e "$(date)\ncram file $D2/$D1 found\n" &>> $CWD/$SM/log/$SM.run.log
        else
            echo -e "$(date)\ncram file $D2/$D1 not found\n", exiting &>> $CWD/$SM/log/$SM.run.log
            scancel -n $SM
        fi
    else
        if [[ -e $D1/$R1.gz && -e $D2/$R2.gz ]]; then
            echo -e "$(date)\nfastq files $D1/$R1.gz and $D2/$R2.gz found\n" &>> $CWD/$SM/log/$SM.run.log
        else
            echo -e "$(date)\nfastq files $D1/$R1.gz and $D2/$R2.gz found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
            scancel -n $SM
        fi
    fi
done


# check if ref exists
echo -e "\n\n$(date)\nChecks for reference files\n" &>> $CWD/$SM/log/$SM.run.log
if [[ -e $REF ]]; then
    echo -e "$(date)\nref file $REF found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nref file $REF not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

if [[ -e $REF.amb && -e $REF.ann && -e $REF.bwt && -e $REF.pac && -e $REF.sa ]]; then
    echo -e "$(date)\nbwa index files found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nbwa index files not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

if [[ -e $REF.dict ]]; then
    echo -e "$(date)\ndict file found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\ndict file not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

if [[ -e $REF.fai ]]; then
    echo -e "$(date)\nfai file found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nfai file not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

# check if targets exist
if [[ -d ${REF%/*}/target_loci/ ]]; then
    echo -e "$(date)\ntarget_loci dir found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\ntarget_loci dir not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi 

# check if recal db exists if recal is required
echo -e "\n\n$(date)\nChecks for BQSR\n" &>> $CWD/$SM/log/$SM.run.log
if [[ $RECAL = true ]]; then
    if [[ -e $RECAL ]]; then
        echo -e "$(date)\nBQSR DB found\n" &>> $CWD/$SM/log/$SM.run.log
    else
        echo -e "$(date)\nBQSR DB required and not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
        scancel -n $SM
    fi
elif [[ $RECAL = false ]]; then
    echo -e "$(date)\nBQSR is not required from user\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(date)\nBQSR var is $REACL, requires true/false, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

# check for final destinations and create some tests
# make the dirs for final destinations 
echo -e "\n\n$(date)\nChecks for final destinations\n" &>> $CWD/$SM/log/$SM.run.log
if [[ -d $BAM ]]; then
    echo -e "$(data)\nFinal bam dir, $BAM, found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(data)\nFinal bam dir, $BAM, not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

if [[ -d $GVCF ]]; then
    echo -e "$(data)\nFinal GVCF dir, $GVCF, found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(data)\nFinal GVCF dir, $GVCF, not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

if [[ -d $LOG ]]; then
    echo -e "$(data)\nFinal log dir, $LOG, found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(data)\nFinal log dir, $LOG, not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

if [[ -d $METRICS ]]; then
    echo -e "$(data)\nFinal metrics dir, $METRICS, found\n" &>> $CWD/$SM/log/$SM.run.log
else
    echo -e "$(data)\nFinal metrics dir, $METRICS, not found, exiting\n" &>> $CWD/$SM/log/$SM.run.log
    scancel -n $SM
fi

module load $BWAMOD
module load $FASTQCMOD 
module load $PIGZMOD
module load $JAVAMOD

java -version
echo $?

fastqc -version
echo $?

pigz -V
echo $?

bwa &> BWAfile.version; grep Version BWAfile.version &>> $CWD/$SM/log/$SM.run.log ; rm BWAfile.version
E=$(bwa | echo $?) # seems to give back a one
# E will carry the exit status

# what happens if a program does not work:
# version will return nothing or an error
# what happens if we call a program that does not exist

echo $? # gets the exit status of the last command. This way we can test if there is a failure.


java -jar $GATK -version
echo $?


java -jar $PICARD MarkDuplicates --versi &> jav.v
# piccard don't work so well
# need to check if result is one line long


$SAMTOOLS --version ; echo $?


if [[ $(echo $?) > 0 ]]; then echo yes; else echo no ; fi





echo -e "\n\n$(date)\nFile and program checks complete, moving on...\n\n\n\n" &>> $CWD/$SM/log/$SM.run.log

