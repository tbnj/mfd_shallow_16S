#!/usr/bin/env bash
#set up environment
set -eu
module purge
module load ea-utils parallel
module load fastp/0.23.2-GCC-10.2.0

mkdir 230707_A00595_0258_BHWJWWDSX5
cd 230707_A00595_0258_BHWJWWDSX5

#make directories
mkdir sequences_trim
mkdir log
mkdir log/json
mkdir log/html

#paths
dir_in=('/incoming/microflora_danica/basecalled/230707_A00595_0258_BHWJWWDSX5')
log_json=('log/json')
log_html=('log/html')
dir_out=('sequences_trim')
tmp=('/projects/microflora_danica/sub_projects/phylotables/tmp')

# find directories
find $dir_in -mindepth 2 -maxdepth 2 -type d ! -path "$dir_in/Stats/*" ! -path "$dir_in/Reports/*" > directories.txt

#rename and create commands for parallel
while read -r line; do
  #echo $line
  dir_name=$(echo ${line##*/} | sed 's/_/-/g' | sed -E 's/[0-9]{5}//')
  #echo $dir_name
  input_R1=$(echo $line/*_R1_001.fastq.gz)
  #echo $input_R1
  input_R2=$(echo $line/*_R2_001.fastq.gz)
  #echo $input_R2
  new_name=$(echo ${input_R1##*/} | sed -E 's/_[^L]+L0/_/' | sed -E 's/_R.+//' | sed -E "s/sample/$dir_name/")
  #echo $new_name
  echo fastp \
  --in1 $input_R1 \
  --in2 $input_R2 \
  --out1 $dir_out/$new_name'_R1.fastq' \
  --out2 $dir_out/$new_name'_R2.fastq' \
  --correction \
  --detect_adapter_for_pe \
  --cut_right \
  --cut_right_window_size 4 \
  --cut_right_mean_quality 20 \
  --average_qual 30 \
  --length_required 100 \
  --dedup \
  --dup_calc_accuracy 6 \
  --thread 6 \
  --overrepresentation_analysis \
  --json $log_json/$new_name'.json' \
  --html $log_html/$new_name'.html' \
  >> command.txt
done < directories.txt

# run commands
cat command.txt | parallel -j5 --tmpdir $tmp &> >(tee -a trim.log)

# multithreaded zipping with pigz
module purge
module load pigz/2.4-foss-2018a

for i in $dir_out/*.fastq; do
  pigz -9 -p 30 $i
done
