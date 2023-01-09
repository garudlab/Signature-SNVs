#!/bin/bash

# ./Step1_qsub_get_file_lengths.sh SuezMontassier_reference_info/species_list.txt 180 SuezMontassier


first_count_input=1
COUNTER="$(($first_count_input-1))"
while IFS= read -r strain;
do
	COUNTER=$((COUNTER + 1)); 
  	echo "$strain" > inputs/data_$COUNTER.in; 
done < $1


qsub -cwd -V -N line_count -l h_data=2G,time=24:00:00 -b y -t $first_count_input:$COUNTER "./Step0_run_get_file_lengths.sh"



