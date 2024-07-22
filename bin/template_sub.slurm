#!/bin/bash
#SBATCH --job-name=JOBNAME              # Job name
#SBATCH --partition=nodes               # What partition the job should run on
#SBATCH --time=TIME                     # Time limit (DD-HH:MM:SS)
#SBATCH --ntasks=NbTASKS                # Number of MPI tasks to request
#SBATCH --cpus-per-task=NPN             # Number of CPU cores per MPI task
#SBATCH --mem=2G                        # Total memory to request
#SBATCH --account=pet-interplan-2024    # Project account to use
#SBATCH --mail-type=END,FAIL            # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=cpm531@york.ac.uk   # Where to send mail
#SBATCH --output=%x-%j.log              # Standard output log
#SBATCH --error=%x-%j.err               # Standard error log

##set -x

TEMPDIR="$HOME/LatHyS/bin/ephemeral/JOBNAME"

SCRATCHDIR="/users/cpm531/scratch/"

cd $SCRATCHDIR
mkdir JOBNAME
cd JOBNAME
mkdir Code
cp -rf $TEMPDIR/src Code/
mkdir Code/Soumission
cp $TEMPDIR/sub_Lathys.slurm Code/Soumission
cp $TEMPDIR/launcher.sh Code/Soumission
cp $TEMPDIR/sub_restart.slurm Code/Soumission
cp ~/LatHyS/bin/restart.sh Code/Soumission

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
mv $HOME/LatHyS/bin/JOBNAME*.err ./
mv $HOME/LatHyS/bin/JOBNAME*.log ./

### Once everything is done, delete the tempdir ###
cd $TEMPDIR
cd ..
rm -rf JOBNAME

### And launch a run on a single node to post-process run results ###
cd $SCRATCHDIR/JOBNAME/ncfiles/
cp $HOME/LatHyS/bin/diagnostic_files/* ./

#source ./launch_diag.sh



