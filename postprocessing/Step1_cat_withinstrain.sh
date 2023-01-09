#!/bin/bash
# 
COUNTER=0


while read -r fam strain;
do
	
    echo $strain;

    signature_snvs/Signature_SNVs_"$strain"*$fam*counts.bz2  > signature_snvs_merged/Signature_SNVs_"$strain"_$fam_counts.bz2 

done < inputs/data_$SGE_TASK_ID.in

