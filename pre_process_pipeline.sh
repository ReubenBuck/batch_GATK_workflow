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
RMOD=R/R-3.3.3

# need to read in the config file! 
while [[ "$1" =~ ^- && ! "$1" == "--" ]]; do case $1 in
-h | --help )
echo "check everything is ok" 
exit
;;
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
runLen=$(expr $(wc -l $CONFIG | cut -d" " -f1) - 1)
;;
-r | --recal )
BQSR=true
;;
-P | --perform )
PERFORM=true
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
-r | --rversion )
shift; RMOD=$1
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


# prepare dirs
sbatch \
--mem=10g --time=2-00:00 --nodes=1 --ntasks=1 \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL,BEGIN --output=prepare_dirs-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/prepare_dirs.sh --sample $SM \
--platform $PL --flowcell $FC --lane $LN --library $LB \
--read1 $R1 --read2 $R2 --path1 $D1 --path2 $D2 \
--ref $REF --workdir $CWD --recal $RECAL --bam $BAM \
--gvcf $GVCF --metrics $METRICS --log $LOG \
--bwa $BWAMOD --fastqc $FASTQCMOD --pigz $PIGZMOD --java $JAVAMOD \
--samtools $SAMTOOLSMOD --gatk $GATK --picard $PICARD \
--runLen $runLen --perform $PERFORM --bqsr $BQSR \

# prepare reads
# this may use some samtools options
sbatch \
--mem=10g --time=2-00:00 --nodes=1 --ntasks=10 --array=1-$runLen \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/prepare_reads-${SM}-%A-%a.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/prepare_reads.sh --sample $SM \
--read1 $R1 --read2 $R2 --path1 $D1 --path2 $D2 --threads 25 \
--workdir $CWD --fastqc $FASTQCMOD --pigz $PIGZMOD \
--samtools $SAMTOOLSMOD --perform $PERFORM \

# Map reads
# bwa does not offer mem usage options but scales proportinally with threads
sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 --array=1-$runLen \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/map_reads-${SM}-${TASKS}-%A-%a.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/map_reads.sh --sample $SM \
--read1 $R1 --read2 $R2 --threads 25 \
--workdir $CWD --bwa $BWAMOD --samtools $SAMTOOLSMOD --perform $PERFORM \
--flowcell $FC --lane $LN --library $LB --platform $PL --ref $REF \

# merge_sort_bams
# for samtools, max mem per task can also be set
# need to go mem/ntask to calculate
# all calculation should be in G
sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/merge_sort_bams-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/merge_sort_bams.sh --sample $SM \
--threads 25 --runLen $runLen \
--workdir $CWD --samtools $SAMTOOLSMOD --perform $PERFORM \

# mark_duplicates.sh
sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/mark_duplicates-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/mark_duplicates.sh --sample $SM \
--workdir $CWD --picard $PICARD --java $JAVAMOD --perform $PERFORM \

# index.sh
IDXJOB=$(sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/index-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/index.sh --sample $SM \
--workdir $CWD --samtools $SAMTOOLSMOD --perform $PERFORM --threads 25 | cut -f 4 -d ' ')

echo $IDXJOB

# unmapped_reads.sh
sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 \
--job-name=${SM}-unmapped --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d afterok:$IDXJOB \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/unmapped_reads-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/unmapped_reads.sh --sample $SM \
--workdir $CWD --samtools $SAMTOOLSMOD --perform $PERFORM --threads 25 \

# realigner_target_creator
# need to check nt and nct defs to determine how to set max mem usage for java 
sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/realigner_target_creator-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/realigner_target_creator.sh --sample $SM \
--workdir $CWD --gatk $GATK --java $JAVAMOD --ref $REF --perform $PERFORM --threads 25 \

# indel realignment 
lociLen=$(ls ${REF%/*}/target_loci | wc -l)
LOCI=$(echo $(ls ${REF%/*}/target_loci) |  sed 's/ /,/g')

# needs a job id for cat_sort_index_bams.sh
CATBAMID=$(sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=2 --array=1-$lociLen \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/indel_realigner-${TASKS}-%A-%a.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/indel_realigner.sh --sample $SM \
--workdir $CWD --gatk $GATK --java $JAVAMOD --ref $REF --perform $PERFORM --loci $LOCI | cut -f 4 -d ' ')

#cat_sort_index.sh
#this is the end bam for non recalibration
#we may not need this at all
#or maybe we can use it later



if [[ $BQSR = true ]]; then
CATBAMID=""

#first_pass_bqsr.s
BQSRID=$(sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/first_pass_bqsr-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/first_pass_bqsr.sh --sample $SM \
--workdir $CWD --gatk $GATK --java $JAVAMOD --recal $RECAL --ref $REF \
--perform $PERFORM --threads 25 | cut -f 4 -d ' ')

echo $BQSRID

# print_reads.sh
# needs a job id for cat_sort_index_bams.sh
CATBAMID=$(sbatch \
--mem=10g --time=2-00:00 --nodes=1 --ntasks=2 --array=1-$lociLen \
--job-name=$SM --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d singleton \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/print_reads-${SM}-%A-%a.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/print_reads.sh --sample $SM \
--workdir $CWD --gatk $GATK --java $JAVAMOD --ref $REF --perform $PERFORM --loci $LOCI | cut -f 4 -d ' ')


# second_pass_bqsr.sh
SECONDBQSRID=$(sbatch \
--mem=50g --time=2-00:00 --nodes=1 --ntasks=10 -d afterok:$BQSRID \
--job-name=${SM}-recal-plots --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/second_pass_bqsr-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/second_pass_bqsr.sh --sample $SM \
--workdir $CWD --gatk $GATK --java $JAVAMOD --recal $RECAL --ref $REF \
--perform $PERFORM --threads 10 --rversion $RMOD | cut -f 4 -d ' ')

echo $SECONFBQSRID

fi

echo $CATBAMID

# cat_sort_index.sh
# this can run parallele to haplotype caller
# needs to take two alternative jobIDs, no we can just overwrite the variable if bqsr is true
# push the output to final destination
# then get an md5sum
sbatch \
--mem=100g --time=2-00:00 --nodes=1 --ntasks=10 \
--job-name=${SM}-cat-bams --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE -d afterok:$CATBAMID \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/cat_sort_index_bams-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/cat_sort_index_bams.sh --sample $SM \
--ref $REF --workdir $CWD --samtools $SAMTOOLSMOD --perform $PERFORM --threads 25 --bqsr $BQSR \


# haplotype caller
haplotypecallerNTASKS=20
haplotypecallerMEM="200g"
haplotypecallerTIME="2-00:00"
sbatch \
--mem=$haplotypecallerMEM --time=$haplotypecallerTIME --ntasks=$haplotypecallerNTASKS \
-d singleton --array 1-$lociLen --nodes=1 \
--job-name=${SM} --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/haplotypecaller-${SM}-%A-%a.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/haplotypecaller.sh --sample $SM \
--workdir $CWD --gatk $GATK --java $JAVAMOD --ref $REF --loci $LOCI --bqsr $BQSR \
--perform $PERFORM --threads $haplotypecallerNTASKS

# cat gvcfs
# for samtools tasks we could put more into mem
cat_gvcfNTASKS=10
cat_gvcfMEM="100g"
cat_gvcfTIME="2-00:00"
sbatch \
--mem=${cat_gvcfMEM} --time=${cat_gvcfTIME} --ntasks=${cat_gvcfNTASKS} \
-d singleton --nodes=1 --job-name=${SM} \
--account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/cat_gvcf-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/cat_gvcf.sh --sample $SM \
--workdir $CWD --picard $PICARD --java $JAVAMOD --ref $REF --loci $LOCI \
--perform $PERFORM --threads ${cat_gvcfNTASKS}


# move everything to final place
FINALID=$(sbatch \
--mem=10g --time=2-00:00 --nodes=1 --ntasks=10 -d singleton \
--job-name=${SM} --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL --output=$CWD/$SM/log/cp_files-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/cp_files.sh --sample $SM \
--workdir $CWD --bam $BAM --metrics $METRICS --log $LOG --gvcf $GVCF --bqsr $BQSR | cut -f 4 -d ' ')

# clean the working dir
if [[ $BQSR = true ]]; then

sbatch \
--mem=10g --time=2-00:00 --nodes=1 --ntasks=1 -d afterok:${FINALID}:${SECONDBQSRID} \
--job-name=${SM} --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL,END --output=clean_wd-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/clean_wd.sh --sample $SM \
--workdir $CWD

else

sbatch \
--mem=10g --time=2-00:00 --nodes=1 --ntasks=1 -d afterok:${FINALID} \
--job-name=${SM} --account=$ACCOUNT --partition=$PARTITION $EXCLUSIVE \
--mail-user=$EMAIL --mail-type=FAIL,END --output=clean_wd-${SM}-%j.out \
/home/buckleyrm/scripts/batch_GATK_workflow/tasks/clean_wd.sh --sample $SM \
--workdir $CWD

fi


