#!/bin/bash
#SBATCH -J JOBNAME     # Request name
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=NODES
#SBATCH --ntasks=NbTASKS
#SBATCH --threads-per-core=1

#SBATCH --constraint=HSW24
#SBATCH --time=TIME 
#SBATCH --output JOBNAME.output

#SBATCH --mail-type=END
#SBATCH --mail-user=clement.moissard@lpp.polytechnique.fr

##set -x


cd ${SCRATCHDIR}
cd JOBNAME
cd ncfiles
mv restart_files/*.nc .

INPUT=$(ls | sort | head -1)
dd=$(echo $INPUT| cut -d'_' -f 1)
mm=$(echo $INPUT| cut -d'_' -f 2)
yy=$(echo $INPUT| cut -d'_' -f 3 | tr -dc '0-9')
date=$(echo $dd'_'$mm'_'$yy)

### on fait le calcul proprement dit
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./quiet_plasma -ncxyz NX NY NZ -dt DT --gstep DX -nhm NHM -pn PLANETNAME -r 1 -rf RF -d $date -ug 1

###Gather all the restart files
mv r_*.nc restart_files

### return in the job root directory ###
cd ..

### Save the .output file separately so that it is easier to access
mv $HOME/bin/JOBNAME.output ./

### Once everything is done, delete the tempdir ###
cd $TEMPDIR
cd ..
rm -rf JOBNAME

### And launch a run on a single node to post-process run results ###
cd ${SCRATCHDIR}/JOBNAME/ncfiles/
cp $HOME/bin/diagnostic_files/* ./

# source ./launch_diag.sh


