#!/bin/bash

### Short description of the job ###

echo "can we make venus?"

JOBNAME='24_08_27_steady_venus_0'

NbTASKS=128                   # Number of tasks to use (MPI processes)

TIME="0-03:15:00"

NX=100
NY=100
NZ=100

TMAX=100
DT=0.02
DX=1

MEM=100

NHM=$(echo $TMAX \/ $DT |bc)

PLANETNAME='venus'
WITH_PLANET='yes'

RESTART=0

### Setup the actual submission file
var1=s/JOBNAME/$JOBNAME/g
var2=s/NODES/$NODES/g
var3=s/MEM/$MEM/g
var4=s/NbTASKS/$NbTASKS/g
var6=s/TIME/$TIME/g
var7=s/NX/$NX/g
var8=s/NY/$NY/g
var9=s/NZ/$NZ/g
var10=s/DT/$DT/g
var11=s/DX/$DX/g
var12=s/NHM/$NHM/g
var13=s/PLANETNAME/$PLANETNAME/g
var14=s/RESTART/$RESTART/g

### Create the sub file, i.e. the file used by the cluster to run the simulation
sed -e ${var1} -e ${var2} -e ${var3} -e ${var4} -e ${var6} -e ${var7} -e ${var8} -e ${var9} -e ${var10} -e ${var11} -e ${var12} -e ${var13} -e ${var14} <template_sub.slurm >sub_Lathys.slurm

### Create the sub_restart file now, so if you need it, you won't have to create### it manually
sed -e ${var1} -e ${var2} -e ${var4} -e ${var6} -e ${var7} -e ${var8} -e ${var9} -e ${var10} -e ${var11} -e ${var12} -e ${var13} -e ${var14} <template_restart.slurm >sub_restart.slurm

### Make def_tregister.F90 start at the beginning
cd ~/LatHyS/Lathys/src/
LASTDUMP=10
var1=s/LASTDUMP/$LASTDUMP/g
sed -e ${var1} <./defs_tregister_modify_this_one.F90 >./defs_tregister.F90

### Compile the code
if [ $WITH_PLANET = 'yes' ]; then
    var=s/YoNPLANET/'#'/g
else
    var=s/YoNPLANET/''/g
fi

sed -e ${var} <~/LatHyS/Lathys/sav_Makefile >~/LatHyS/Lathys/Makefile

cd ~/LatHyS/Lathys/
make clean
make
make diag
cd ~/LatHyS/bin/

### Moving the files before entering the queue, so that I can launch ###
### multiple runs at the same time ###
TEMPDIR="./ephemeral/$JOBNAME"
if [ -d $TEMPDIR ]
then
	rm -rf $TEMPDIR
fi 
mkdir $TEMPDIR
## Bring the necessary files to the tempdir ##
cp -rf $HOME/LatHyS/Lathys/src $TEMPDIR/
cp -rf $HOME/LatHyS/Lathys/quiet_plasma $TEMPDIR/
cp -rf $HOME/LatHyS/Lathys/diag $TEMPDIR/
mv $HOME/LatHyS/bin/sub_Lathys.slurm $TEMPDIR/
cp $HOME/LatHyS/bin/launcher.sh $TEMPDIR/
mv $HOME/LatHyS/bin/sub_restart.slurm $TEMPDIR/

sbatch $TEMPDIR/sub_Lathys.slurm



