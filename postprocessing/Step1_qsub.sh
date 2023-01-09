#!/bin/bash

# ./Step1_qsub.sh species_list.txt family_list.txt


COUNTER=1
mkdir signature_snvs_merged

while IFS= read -r strain;
do
	while IFS= read -r fam;
	do
		COUNTER=$((COUNTER + 1));        
		echo $strain;
		echo $fam
		echo "$fam $strain" > inputs/data_$COUNTER.in; 
	done < $2
done < $1

qsub -cwd -V -N MergeSignatureSNV -e misc -o misc -l h_data=6G,time=10:00:00 -b y -t $first_count_input:$COUNTER "./Step1_cat_withinstrain.sh"
