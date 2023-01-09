#!/bin/bash

mkdir signature_snvs_per_family

while read -r fam ;
do

    cat signature_snvs_merged/Signature_SNVs*$fam_counts.csv.bz2 >> signature_snvs_per_family/Signature_SNVs_$fam_counts.csv.bz2
    

done < inputs/data_$SGE_TASK_ID.in

