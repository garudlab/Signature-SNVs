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

outputnum1=$1
outputnum2=$2
# echo job info on joblog:
echo "Name $outputnum1"


# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load anaconda3/2020.11

# activate an already existing conda environment (CHANGE THE NAME OF THE ENVIRONMENT):
conda activate /u/home/b/briscoel/project-ngarud/inStrain

inStrain compare -i ../Mom"$outputnum1"/mom_profile_"$outputnum1" ../Baby"$outputnum2"/baby_profile_"$outputnum2"