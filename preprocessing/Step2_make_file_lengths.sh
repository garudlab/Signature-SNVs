#!/bin/bash

# ./Step2_make_file_lengths.sh species_list.txt

while IFS= read -r strain;
do
	COUNTER=$((COUNTER + 1)); 
	var=$( cat species_file_lengths/"$strain"_lengths.txt)
  	echo "$var" >> species_lengths.txt; 
done < $1
