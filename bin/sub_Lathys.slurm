#!/bin/bash
#SBATCH -J 23_01_19_Long_box_run_0     # Request name
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --threads-per-core=1
#SBATCH --dependency=singleton
#SBATCH --constraint=HSW24
#SBATCH --time=23:30:00 
#SBATCH --output 23_01_19_Long_box_run_0.output
#SBATCH --mail-type=END
#SBATCH --mail-user=clement.moissard@lpp.polytechnique.fr

##set -x

TEMPDIR="$HOME/bin/ephemeral/23_01_19_Long_box_run_0"

cd ${SCRATCHDIR}
mkdir 23_01_19_Long_box_run_0
cd 23_01_19_Long_box_run_0
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
srun --mpi=pmi2 -K1 --resv-ports -n $SLURM_NTASKS ./quiet_plasma -ncxyz 1000 8 8 -dt 0.05 --gstep 1 -nhm 6000 -pn earth -r 0 -ug 1

###Gather all the restart files
mkdir restart_files
mv r_*.nc restart_files

### return in the job root directory ###
cd ..

### Save the .output file separately so that it is easier to access
mv $HOME/bin/23_01_19_Long_box_run_0.output ./
mv $HOME/bin/23_01_19_Long_box_run_0_post.output ./

### Once everything is done, delete the tempdir ###
cd $TEMPDIR
cd ..
rm -rf 23_01_19_Long_box_run_0

### And launch a run on a single node to post-process run results ###
cd ${SCRATCHDIR}/23_01_19_Long_box_run_0/ncfiles/
cp $HOME/bin/diagnostic_files/* ./

#source ./launch_diag.sh



