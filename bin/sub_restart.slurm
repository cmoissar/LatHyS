#!/bin/bash
#SBATCH -J 23_01_19_Long_box_run_0     # Request name
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --threads-per-core=1

#SBATCH --constraint=HSW24
#SBATCH --time=23:30:00 
#SBATCH --output 23_01_19_Long_box_run_0.output

#SBATCH --mail-type=END
#SBATCH --mail-user=clement.moissard@lpp.polytechnique.fr

##set -x


cd ${SCRATCHDIR}
cd 23_01_19_Long_box_run_0
cd ncfiles
mv restart_files/*.nc .

INPUT=$(ls | sort | head -1)
dd=$(echo $INPUT| cut -d'_' -f 1)
mm=$(echo $INPUT| cut -d'_' -f 2)
yy=$(echo $INPUT| cut -d'_' -f 3 | tr -dc '0-9')
date=$(echo $dd'_'$mm'_'$yy)

### on fait le calcul proprement dit
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./quiet_plasma -ncxyz 1000 8 8 -dt 0.05 --gstep 1 -nhm 6000 -pn earth -r 1 -rf RF -d $date -ug 1

###Gather all the restart files
mv r_*.nc restart_files

### return in the job root directory ###
cd ..

### Save the .output file separately so that it is easier to access
mv $HOME/bin/23_01_19_Long_box_run_0.output ./

### Once everything is done, delete the tempdir ###
cd $TEMPDIR
cd ..
rm -rf 23_01_19_Long_box_run_0

### And launch a run on a single node to post-process run results ###
cd ${SCRATCHDIR}/23_01_19_Long_box_run_0/ncfiles/
cp $HOME/bin/diagnostic_files/* ./

# source ./launch_diag.sh


