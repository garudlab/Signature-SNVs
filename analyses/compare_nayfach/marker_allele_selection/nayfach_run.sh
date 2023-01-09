#!/bin/bash
#$ -N Nayfach
#$ -e logs
#$ -o logs
#$ -cwd
#$ -r y
#$ -j y
#$ -b y
#$ -l h_data=12G
#$ -l time=23:00:00
#$ -tc 100 # Throttle to max 100 tasks at a time
#$ -t 1-95

 
. /u/local/Modules/default/init/modules.sh
module load anaconda3/2020.11


readarray species < /u/home/b/briscoel/project-ngarud/FEASTX/SimulationDFiles/species_lists/species_list.txt
species=(null ${species[@]}) # zero to one start index
specie=${species[$SGE_TASK_ID]}
echo $specie      
python nayfach_main.py --strain $specie --seed 27

