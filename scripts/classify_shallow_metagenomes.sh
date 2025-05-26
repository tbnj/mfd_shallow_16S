#!/usr/bin/env bash
set -eu
module purge

#### Classification of 16S and 18S fragments from metagenomes     ####
#### The pipeline relies on a functioning installation of USEARCH ####
#### Run pipeline for each individual sequencing run (RUNNAME)    ####

## Function to add a header to echo, for a better console output overview
echoWithHeader() {
  echo " *** [$(date '+%Y-%m-%d %H:%M:%S')]: $*"
}

## Function to report total run time 
echoDuration() {
  duration=$(printf '%02dh:%02dm:%02ds\n' $((SECONDS/3600)) $((SECONDS%3600/60)) $((SECONDS%60)))
  echoWithHeader "Done in: $duration!"
}

### Set up environment
DIR_IN=('/mfd_shallow_16S/data/RUNNAME/')
SEQ_IN=('/mfd_shallow_16S/data/RUNNAME/sequences_trim/')
DIR_OUT_F=('forward')
DIR_OUT_R=('reverse')
HMMS=('/mfd_shallow_16S/databases/HMMS')
UDB=('/mfd_shallow_16S/databases/MFG_ssu_database_v1.3_NR987_sintax.udb')
THREADS_USEARCH=60
THREADS_NHMMER=echo "$((THREADS_USEARCH/3))"

## Create directories
cd $DIR_IN
mkdir sintax_classification
cd sintax_classification

## List of sample names for classification
find $SEQ_IN/*_R1.fastq.gz > samples_classification.txt

## Run pipleine on list of samples (handles both R1 and R2)
while read -r line; do
  INPUT=$(echo ${line} | sed -E 's/_R.+/_/')
  NEW_DIR=$(echo ${line##*/} | sed -E 's/_R.*//')
  NEW_NAME=$(echo ${line##*/} | sed -E 's/_[^-]+$//')
  echoWithHeader "  - Starting analysis of sample $name..."
  mkdir $NEW_DIR $NEW_DIR/forward $NEW_DIR/reverse $NEW_DIR/tmp
  echoWithHeader "  - Searching forward reads for 16S and 18S fragments..."
  module purge
  module load HMMER/3.3.2-foss-2020b
  zcat $INPUT'R1.fastq.gz' | \
  awk '{print ">" substr($0,2);getline;print;getline;getline}' - | \
  tee >(nhmmer --incE 1e-05 -E 1e-05 --cpu THREADS_NHMMER -o /dev/null --noali --tblout $NEW_DIR/forward/'bac_'$NEW_NAME'_forward.hmmout.txt' $HMMS/bac.hmm -) | \
  tee >(nhmmer --incE 1e-05 -E 1e-05 --cpu THREADS_NHMMER -o /dev/null --noali --tblout $NEW_DIR/forward/'arc_'$NEW_NAME'_forward.hmmout.txt' $HMMS/arc.hmm -) | \
  nhmmer --incE 1e-05 -E 1e-05 --cpu THREADS_NHMMER -o /dev/null --noali --tblout $NEW_DIR/forward/'euk_'$NEW_NAME'_forward.hmmout.txt' $HMMS/euk.hmm -
  echoWithHeader "  - Searching reverse reads for 16S and 18S fragments..."
  zcat $INPUT'R2.fastq.gz' | \
  awk '{print ">" substr($0,2);getline;print;getline;getline}' - | \
  tee >(nhmmer --incE 1e-05 -E 1e-05 --cpu THREADS_NHMMER -o /dev/null --noali --tblout $NEW_DIR/reverse/'bac_'$NEW_NAME'_reverse.hmmout.txt' $HMMS/bac.hmm -) | \
  tee >(nhmmer --incE 1e-05 -E 1e-05 --cpu THREADS_NHMMER -o /dev/null --noali --tblout $NEW_DIR/reverse/'arc_'$NEW_NAME'_reverse.hmmout.txt' $HMMS/arc.hmm -) | \
  nhmmer --incE 1e-05 -E 1e-05 --cpu THREADS_NHMMER -o /dev/null --noali --tblout $NEW_DIR/reverse/'euk_'$NEW_NAME'_reverse.hmmout.txt' $HMMS/euk.hmm -
  echoWithHeader "  - Extracting reads"
  module purge 
  awk -F " " 'NR>2 {print $1}' $NEW_DIR/forward/'bac_'$NEW_NAME'_forward.hmmout.txt' | grep -vE "^#" > $NEW_DIR/tmp/forward_IDs.txt
  awk -F " " 'NR>2 {print $1}' $NEW_DIR/forward/'arc_'$NEW_NAME'_forward.hmmout.txt' | grep -vE "^#" >> $NEW_DIR/tmp/forward_IDs.txt
  awk -F " " 'NR>2 {print $1}' $NEW_DIR/forward/'euk_'$NEW_NAME'_forward.hmmout.txt' | grep -vE "^#" > $NEW_DIR/tmp/forward_euk_IDs.txt
  sort $NEW_DIR/tmp/forward_IDs.txt | uniq - > $NEW_DIR/tmp/forward_IDs_unique.txt
  awk -F " " 'NR>2 {print $1}' $NEW_DIR/reverse/'bac_'$NEW_NAME'_reverse.hmmout.txt' | grep -vE "^#" > $NEW_DIR/tmp/reverse_IDs.txt
  awk -F " " 'NR>2 {print $1}' $NEW_DIR/reverse/'arc_'$NEW_NAME'_reverse.hmmout.txt' | grep -vE "^#" >> $NEW_DIR/tmp/reverse_IDs.txt
  awk -F " " 'NR>2 {print $1}' $NEW_DIR/reverse/'euk_'$NEW_NAME'_reverse.hmmout.txt' | grep -vE "^#" > $NEW_DIR/tmp/reverse_euk_IDs.txt
  sort $NEW_DIR/tmp/reverse_IDs.txt | uniq - > $NEW_DIR/tmp/reverse_IDs_unique.txt
  zcat $INPUT'R1.fastq.gz' | grep -A 3 -f --no-group-separator $NEW_DIR/tmp/forward_IDs_unique.txt > $NEW_DIR/forward/'arc_bac_'$NEW_NAME'_forward.fq'
  zcat $INPUT'R1.fastq.gz' | grep -A 3 -f --no-group-separator $NEW_DIR/tmp/forward_euk_IDs.txt > $NEW_DIR/forward/'euk_'$NEW_NAME'_forward.fq'
  zcat $INPUT'R2.fastq.gz' | grep -A 3 -f --no-group-separator $NEW_DIR/tmp/forward_IDs_unique.txt > $NEW_DIR/reverse/'arc_bac_'$NEW_NAME'_reverse.fq'
  zcat $INPUT'R2.fastq.gz' | grep -A 3 -f --no-group-separator $NEW_DIR/tmp/forward_euk_IDs.txt > $NEW_DIR/reverse/'euk_'$NEW_NAME'_reverse.fq'
  echoWithHeader "  - Classifying reads"
  module purge
  usearch -sintax $NEW_DIR/forward/'arc_bac_'$NEW_NAME'_forward.fq' -db $UDB -tabbedout $NEW_DIR/forward/'arc_bac_'$NEW_NAME'_forward_MFG_ssu_database_NR987_trunc.sintax' \
  -strand both -sintax_cutoff 0.8 -threads THREADS_USEARCH -quiet
  usearch -sintax $NEW_DIR/reverse/'arc_bac_'$NEW_NAME'_reverse.fq' -db $UDB -tabbedout $NEW_DIR/reverse/'arc_bac_'$NEW_NAME'_reverse_MFG_ssu_database_NR987_trunc.sintax' \
  -strand both -sintax_cutoff 0.8 -threads THREADS_USEARCH -quiet
  echoDuration 
done < ../samples_classification.txt &> >(tee -a classification.log)
