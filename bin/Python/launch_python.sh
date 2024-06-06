#!/bin/bash

INPUT=$(pwd)
JOBNAME=$(echo $INPUT| cut -d'/' -f 6)

var=s/JOBNAME/$JOBNAME/g

sed -e ${var} <template_python.slurm >python.slurm

sbatch ./python.slurm
