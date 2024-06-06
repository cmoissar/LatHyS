#!/bin/bash

### Short description of the job ###

echo "Can  I run this code on York's research computers?"

JOBNAME='23_01_19_Long_box_run_0'

NODES=1
NPN=24                        # Keep NPN=24 unless you change --constraint=HSW24 in the template
NbTASKS=$(expr $NODES \* $NPN)   # Number of tasks to use (MPI processes)

#TIME="'2-12:00:00'"
TIME="23:30:00"

NX=1000
NY=8
NZ=8

TMAX=300
DT=0.05
DX=1

NHM=$(echo $TMAX \/ $DT |bc)

PLANETNAME='earth'
WITH_PLANET='no'

RESTART=0

### Setup the actual submission file
var1=s/JOBNAME/$JOBNAME/g
var2=s/NODES/$NODES/g
var3=s/NPN/$NPN/g
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
sed -e ${var1} -e ${var2} -e ${var3} -e ${var4} -e ${var6} -e ${var7} -e ${var8} -e ${var9} -e ${var10} -e ${var11} -e ${var12} -e ${var13} -e ${var14} <template_restart.slurm >sub_restart.slurm

### Make def_tregister.F90 start at the beginning
cd ~/Lathys/src/
LASTDUMP=10
var1=s/LASTDUMP/$LASTDUMP/g
sed -e ${var1} <./defs_tregister_modify_this_one.F90 >./defs_tregister.F90

### Compile the code
if [ $WITH_PLANET = 'yes' ]; then
    var=s/YoNPLANET/'#'/g
else
    var=s/YoNPLANET/''/g
fi

sed -e ${var} <~/Lathys/sav_Makefile >~/Lathys/Makefile

cd ~/Lathys/
make clean
make
cd ~/bin/

### Moving the files before entering the queue, so that I can launch ###
### multiple runs at the same time ###
TEMPDIR="./ephemeral/$JOBNAME"
if [ -d $TEMPDIR ]
then
	rm -rf $TEMPDIR
fi 
mkdir $TEMPDIR
## Bring the necessary files to the tempdir ##
cp -rf $HOME/Lathys/src $TEMPDIR/
cp -rf $HOME/Lathys/quiet_plasma $TEMPDIR/
cp -rf $HOME/Lathys/diag $TEMPDIR/
mv $HOME/bin/sub_Lathys.slurm $TEMPDIR/
cp $HOME/bin/launcher.sh $TEMPDIR/
mv $HOME/bin/sub_restart.slurm $TEMPDIR/

qsub $TEMPDIR/sub_Lathys.slurm



