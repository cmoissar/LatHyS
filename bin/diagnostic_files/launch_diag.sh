#!/bin/bash

mv ./*.txt temporal_B/

INPUT=$(pwd)
JOBNAME=$(echo $INPUT| cut -d'/' -f 6)

var=s/JOBNAME/$JOBNAME/g

sed -e ${var} <template_diag.slurm >diag_Lathys.slurm

sbatch ./diag_Lathys.slurm
