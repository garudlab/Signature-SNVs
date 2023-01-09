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

acc=$1
# echo job info on joblog:
echo "Name $outputname"


# load the job environment:
. /u/local/Modules/default/init/modules.sh
module load anaconda3/2020.11

# activate an already existing conda environment (CHANGE THE NAME OF THE ENVIRONMENT):
conda activate /u/home/b/briscoel/project-ngarud/inStrain

inStrain profile /u/home/b/briscoel/project-ngarud/FEASTX/SimulationDFiles/inStrainOutput/$acc/genomes_"$acc".bam /u/home/b/briscoel/project-ngarud/FEASTX/SimulationDFiles/inStrainOutput/$acc/genomes_"$acc".fa -o /u/home/b/briscoel/project-ngarud/FEASTX/SimulationDFiles/inStrainOutput/$acc/instrain_"$acc"
