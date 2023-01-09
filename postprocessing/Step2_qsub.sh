
#!/bin/bash
#$ -N zipacrossstrains
#$ -e misc
#$ -o misc
#$ -cwd
#$ -l h_data=2G  #20G
#$ -l time=4:00:00
#$ -tc 100 # Throttle to max 100 tasks at a time
#$ -t 1-200 # 200 or number of strains if > 200 
#$ -l highp


readarray families < $1

families=(null ${families[@]}) # zero to one start index
family=${families[$SGE_TASK_ID]}
echo $family

postprocessing/Step2_catzip_acrossstrain.sh $family

