#!/bin/bash
#SBATCH -J diag_big_JOBNAME # Request name
#SBATCH --ntasks-per-node=24
#SBATCH --nodes=1
#SBATCH --ntasks=24
#SBATCH --threads-per-core=1

#SBATCH --constraint=HSW24
#SBATCH --time=20:00:00 
#SBATCH --output Diag.output

#SBATCH --mail-type=END
#SBATCH --mail-user=clement.moissard@lpp.polytechnique.fr

cd ${SCRATCHDIR}/JOBNAME/ncfiles/

source ./script_diag.sh

