#!/bin/bash
#$ -N instrain_job
#$ -e logs
#$ -o logs
#$ -cwd
#$ -r y
#$ -j y
#$ -b y
#$ -l h_data=10G
#$ -l time=48:00:00
#$ -l highp
#$ -pe shared 4
#$ -t 1-280

readarray accs1 < ../accessions_list/instrain_compare1.txt
accs1=(null ${accs1[@]}) # zero to one start index
acc1=${accs1[$SGE_TASK_ID]}
echo $acc1

readarray accs2 < ../accessions_list/instrain_compare2.txt
accs2=(null ${accs2[@]}) # zero to one start index
acc2=${accs2[$SGE_TASK_ID]}
echo $acc2




# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load anaconda3/2020.11

# activate an already existing conda environment (CHANGE THE NAME OF THE ENVIRONMENT):
conda activate /u/home/b/briscoel/project-ngarud/inStrain

inStrain compare -i /u/home/b/briscoel/project-ngarud/FEASTX/SimulationDFiles/inStrainOutput/"$acc1"/instrain_"$acc1" /u/home/b/briscoel/project-ngarud/FEASTX/SimulationDFiles/inStrainOutput/"$acc2"/instrain_"$acc2" -o /u/home/b/briscoel/project-ngarud/FEASTX/SimulationDFiles/inStrainCompare/CompareCV5_"$acc1"_"$acc2" -p 4