#!/bin/bash

### Short description of the job ###

echo "Can I make the launching process simpler?"

JOBNAME='24_05_15_Long_box_run_1'

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
make all
cd ~/bin/

### Moving the files before launching into their own RUN forlder, so that I can launch ###
### multiple runs at the same time ###
RUNDIR="./RUNS/$JOBNAME"
if [ -d $RUNDIR ]
then
	rm -rf $RUNDIR
fi 
mkdir $RUNDIR
## Bring the necessary files to the tempdir ##
cp -rf $HOME/Lathys/src $RUNDIR/
cp -rf $HOME/Lathys/quiet_plasma $RUNDIR/
cp -rf $HOME/Lathys/diag $RUNDIR/
cp $HOME/bin/york_launcher.sh $RUNDIR/

cd $RUNDIR

./quiet_plasma -ncxyz $NX $NY $NZ -dt $DT --gstep $DX -nhm $NHM -pn $PLANETNAME -r $RESTART -ug 1


###Gather all the restart files
mkdir restart_files
mv r_*.nc restart_files

### Copy the files necessary for the diagnostic to the RUN directory ###
cd ncfiles/
cp $HOME/bin/diagnostic_files/script_diag.sh ./

### Concatenate results **
./script_diag.sh

