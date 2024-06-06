#!/bin/bash
#SBATCH -J JOBNAME     # Request name
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=NODES
#SBATCH --ntasks=NbTASKS
#SBATCH --threads-per-core=1
#SBATCH --dependency=singleton
#SBATCH --constraint=HSW24
#SBATCH --time=TIME 
#SBATCH --output JOBNAME.output
#SBATCH --mail-type=END
#SBATCH --mail-user=clement.moissard@lpp.polytechnique.fr

##set -x

TEMPDIR="$HOME/bin/ephemeral/JOBNAME"

cd ${SCRATCHDIR}
mkdir JOBNAME
cd JOBNAME
mkdir Code
cp -rf $TEMPDIR/src Code/
mkdir Code/Soumission
cp $TEMPDIR/sub_Lathys.slurm Code/Soumission
cp $TEMPDIR/launcher.sh Code/Soumission
cp $TEMPDIR/sub_restart.slurm Code/Soumission
cp ~/bin/restart.sh Code/Soumission

mkdir ncfiles
cd ncfiles

#
### on amene sur le workdir les fichiers necessaires
#
cp $TEMPDIR/quiet_plasma .
cp $TEMPDIR/diag .
### on fait le calcul proprement dit
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./quiet_plasma -ncxyz NX NY NZ -dt DT --gstep DX -nhm NHM -pn PLANETNAME -r RESTART -ug 1

###Gather all the restart files
mkdir restart_files
mv r_*.nc restart_files

### return in the job root directory ###
cd ..

### Save the .output file separately so that it is easier to access
mv $HOME/bin/JOBNAME.output ./
mv $HOME/bin/JOBNAME_post.output ./

### Once everything is done, delete the tempdir ###
cd $TEMPDIR
cd ..
rm -rf JOBNAME

### And launch a run on a single node to post-process run results ###
cd ${SCRATCHDIR}/JOBNAME/ncfiles/
cp $HOME/bin/diagnostic_files/* ./

#source ./launch_diag.sh



