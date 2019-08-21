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
--fastqc )
shift; FASTQCMOD=$1
;;
--bwa )
shift; BWAMOD=$1
;;
--java )
shift; JAVAMOD=$1
;;
--perform )
shift; PERFORM=$1
;;
--bqsr )
shift; BQSR=$1
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
    echo -e "$(date)\tbegin\tprepare_dirs.sh" &>> $CWD/$SM/log/$SM.run.log
else
    mkdir -p $CWD/$SM/log
    echo -e "$(date)\tbegin\tprepare_dirs.sh" &>> $CWD/$SM/log/$SM.run.log
fi
mkdir -p $CWD/$SM/fastq
mkdir -p $CWD/$SM/bam
mkdir -p $CWD/$SM/metrics
mkdir -p $CWD/$SM/gvcf
mkdir -p $CWD/$SM/tmp


echo -e "prepare_dirs.sh is running on $(hostname)" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out

# here we can start measuring performance stats
echo -e "\n\n$(date)\nChecks for performance\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
if [[ $PERFORM = true ]]; then
    echo -e "$(date)\nSetting up task performance metrics\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "$(date): prepare_dirs.sh is running on $(hostname)" >  $CWD/$SM/metrics/perform_prepare_dirs_$SM.txt
    vmstat -twn -S m 1 >> $CWD/$SM/metrics/perform_prepare_dirs_$SM.txt &
elif [[ $PERFORM = false ]]; then
    echo -e "$(date)\nPerformance metrics not recorded\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nPerformance var is $PERFORM, requires true/false, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

# Check for programs

echo -e "\n\n$(date)\nCheck for programs\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out

module load $JAVAMOD
module load $SAMTOOLSMOD
module load $FASTQCMOD
module load $PIGZMOD
module load $BWAMOD


java -version; javExit=$?
if [[ $javExit = 0 ]]; then
    echo -e "$(date)\nUsing java version:" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    java -version &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nJava did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

samtools --version; samExit=$?
if [[ $samExit = 0 ]]; then
    echo -e "$(date)\nUsing samtools version:" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    samtools --version &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nSamtools did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

fastqc --version; fasExit=$?
if [[ $fasExit = 0 ]]; then
    echo -e "$(date)\nUsing fastqc version:" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    fastqc --version &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nFastqc did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

pigz -V; pigExit=$?
if [[ $pigExit = 0 ]]; then
    echo -e "$(date)\nUsing pigz version:" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    pigz -V &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\Pigz did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

java -jar $GATK -version; gatExit=$?
if [[ $gatExit = 0 ]]; then
    echo -e "$(date)\nUsing GATK version:" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    java -jar $GATK -version &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\GATK did not exit with 0 status, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

bwa &> $CWD/$SM/tmp/bwa.open.txt
grep Version $CWD/$SM/tmp/bwa.open.txt > $CWD/$SM/tmp/bwa.version.txt
if [[ $(wc -l $CWD/$SM/tmp/bwa.version.txt | cut -d" " -f1) = 1 ]]; then
    echo -e "$(date)\nUsing BWA version:" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    cat $CWD/$SM/tmp/bwa.version.txt &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nBWA did not report version, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi


java -jar $PICARD MarkDuplicates --version &> $CWD/$SM/tmp/picard.version.txt
if [[ $(wc -l $CWD/$SM/tmp/picard.version.txt | cut -d" " -f1) = 1 ]]; then
    echo -e "$(date)\nUsing picard MarkDuplicates version:" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    cat $CWD/$SM/tmp/picard.version.txt &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    echo -e "\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nPicard MarkDuplicates did not report version, exiting" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

echo -e "\n\n$(date)\nProgram checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out


# check files exist!
# check if unaligned files exist
echo -e "\n\n$(date)\nChecks for unaligned WGS data\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
for i in $(seq 1 $runLen); do
    if [[ $D2 = *".bam" ]]; then
        if [[ -e $D1/$D2 ]]; then
            echo -e "$(date)\nbam file $D2/$D1 found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
        else
            echo -e "$(date)\nbam file $D2/$D1 not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
            scancel -n $SM
        fi
    elif [[ $D2 = *".cram" ]]; then
        if [[ -e $D1/$D2 ]]; then
            echo -e "$(date)\ncram file $D2/$D1 found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
        else
            echo -e "$(date)\ncram file $D2/$D1 not found\n", exiting &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
            scancel -n $SM
        fi
    else
        if [[ -e ${D1arr[$i]}/${R1arr[$i]}.gz && -e ${D2arr[$i]}/${R2arr[$i]}.gz ]]; then
            echo -e "$(date)\nfastq files ${D1arr[$i]}/${R1arr[$i]}.gz and ${D2arr[$i]}/${R2arr[$i]}.gz found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
        else
            echo -e "$(date)\nfastq files ${D1arr[$i]}/${R1arr[$i]}.gz and ${D2arr[$i]}/${R2arr[$i]}.gz not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
            scancel -n $SM
        fi
    fi
done

echo -e "\n\n$(date)\nData file checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out


# check if ref exists
echo -e "\n\n$(date)\nChecks for reference files\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
if [[ -e $REF ]]; then
    echo -e "$(date)\nref file $REF found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nref file $REF not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

if [[ -e $REF.amb && -e $REF.ann && -e $REF.bwt && -e $REF.pac && -e $REF.sa ]]; then
    echo -e "$(date)\nbwa index files found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nbwa index files not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

if [[ -e ${REF/.fa/}.dict ]]; then
    echo -e "$(date)\ndict file found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\ndict file not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

if [[ -e $REF.fai ]]; then
    echo -e "$(date)\nfai file found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nfai file not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

# check if targets exist
if [[ -d ${REF%/*}/target_loci/ ]]; then
    echo -e "$(date)\ntarget_loci dir found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\ntarget_loci dir not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi 

# check if recal db exists if recal is required
echo -e "\n\n$(date)\nChecks for BQSR\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
if [[ $BQSR = true ]]; then
    if [[ -e $RECAL ]]; then
        echo -e "$(date)\nBQSR DB found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    else
        echo -e "$(date)\nBQSR DB required and not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
        scancel -n $SM
    fi
elif [[ $BQSR = false ]]; then
    echo -e "$(date)\nBQSR is not required from user\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nBQSR var is $BQSR, requires true/false, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

echo -e "\n\n$(date)\nReference file checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out



# check for final destinations and create some tests
# make the dirs for final destinations 
echo -e "\n\n$(date)\nChecks for final destinations\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
if [[ -d $BAM ]]; then
    echo -e "$(date)\nFinal bam dir, $BAM, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nFinal bam dir, $BAM, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

if [[ -d $GVCF ]]; then
    echo -e "$(date)\nFinal GVCF dir, $GVCF, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nFinal GVCF dir, $GVCF, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

if [[ -d $LOG ]]; then
    echo -e "$(date)\nFinal log dir, $LOG, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nFinal log dir, $LOG, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

if [[ -d $METRICS ]]; then
    echo -e "$(date)\nFinal metrics dir, $METRICS, found\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
else
    echo -e "$(date)\nFinal metrics dir, $METRICS, not found, exiting\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out
    scancel -n $SM
fi

echo -e "\n\n$(date)\nFinal destination checks complete.....\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out

echo -e "\n\n$(date)\nFile and program checks complete, moving on...\n\n\n\n" &>> $CWD/$SM/log/prepare_dirs-${SM}-%j.out

echo -e "$(date)\tend\tprepare_dirs.sh" &>> $CWD/$SM/log/$SM.run.log