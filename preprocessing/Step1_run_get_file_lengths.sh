#!/bin/bash

mkdir species_file_lengths;

while read -r strain; do
	echo $strain;
	path=midas_output/snps/"$strain"/;
	bzcat "$path"snps_ref_freq.txt.bz2 | wc -l > species_file_lengths/"$strain"_lengths.txt;

done < inputs/data_$SGE_TASK_ID.in
