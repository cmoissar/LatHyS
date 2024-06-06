#!/bin/bash
echo "Challenge restart: make it simple"

Run_directory="$(pwd)/../.."

LASTDUMP=???
RF=???  #(0 or 1)
WITH_PLANET=???   #('yes' or 'no')

#TODO: LASTDUMP as keyboard entry
#TODO: streamline the choice of the -rf (0 or 1)

### Make sure the previous output files are not overwritten
output_path=$(ls ../../*.output)
output_name=$(echo $output_path| cut -d'/' -f 3| cut -d'.' -f 1)
mkdir ../../outputfiles
new_name=$(echo $output_name'_a.output')
if test -f $(echo '../../outputfiles/'$new_name); then
    new_name=$(echo $output_name'_b.output')
fi
if test -f $(echo '../../outputfiles/'$new_name); then
    new_name=$(echo $output_name'_c.output')
fi
if test -f $(echo '../../outputfiles/'$new_name); then
    new_name=$(echo $output_name'_d.output')
fi
mv $output_path $(echo '../../outputfiles/'$new_name)
# The 14 last lines are quite ugly.
# The way forward is to use $str = 'a'; echo ++$str; which gives b,
# then c. 

### Modify def_tregister.F90
cd $Run_directory/Code/src/
echo "I am there: $(pwd)"
var1=s/LASTDUMP/$LASTDUMP/g
sed -e ${var1} <./defs_tregister_modify_this_one.F90 >./defs_tregister.F90

### Modify sub_restart.slurm
cd $Run_directory/Code/Soumission/
echo "And now, I am there: $(pwd)"
cp sub_restart.slurm copy_sub_restart.slurm
var2=s/RF/$RF/g
sed -e ${var2} <./copy_sub_restart.slurm >./sub_restart.slurm
echo "Changed def_tregister.F90 and sub_restart.slurm"

### Empty ~/Lathys/src
rm -rf ~/Lathys/src/*
echo "Emptied ~/Lathys/src/ !"

### Replace it with the current run's src
cp -rf $Run_directory/Code/src/* ~/Lathys/src/
echo "And replaced it anew!"

### Compile the code
if [ $WITH_PLANET = 'yes' ]; then
    var=s/YoNPLANET/'#'/g
else
    var=s/YoNPLANET/''/g
fi

sed -e ${var} <~/Lathys/sav_Makefile >~/Lathys/Makefile

cd ~/Lathys
echo "I am at $(pwd) now, and I am going to compile the code. Might take a while, get a coffee!"
make clean
make
make diag
echo "I compiled everything for you, you are welcome"

### Copy the relevant files to the simulation directory
cd $Run_directory/ncfiles
rm -f quiet_plasma
rm -f diag
cp ~/Lathys/quiet_plasma .
cp ~/Lathys/diag .
echo "I am back there: $(pwd)"
echo "And I just replaced quiet_plasma and diag by their new versions"
mv temporal_B/* .
echo "Also made sure that the virtual satellites' data are in the right place"

### Launch the restart
cd $Run_directory/Code/Soumission
echo "I am here: $(pwd)"
rm -f ~/bin/sub_restart.slurm
cp ./sub_restart.slurm ~/bin/
cd ~/bin
qsub ./sub_restart.slurm
echo "I am here: $(pwd)"
echo "and I launched sub_restart.slurm !"

### Reinitialize sub_restart.slurm in case you need it again
cd $Run_directory/Code/Soumission
rm -f sub_restart.slurm
cp -f copy_sub_restart.slurm sub_restart.slurm
rm -f copy_sub_restart.slurm
echo "Back here: $(pwd), again"
echo "Reinitialized sub_restart.slurm"



