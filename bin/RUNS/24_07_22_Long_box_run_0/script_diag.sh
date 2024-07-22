#!/bin/bash

# Script in case something went wrong at the end of the run, but there are still data to play with

###Gather all the restart files
#mkdir restart_files
#mv r_*.nc restart_files

### datas from the simulation are scattered over the procs,
### diag bring them back together
a=300
Dt_diag=1

### Get the date of the simulation, needed for the "diag" command
INPUT=$(ls | sort | head -1)
dd=$(echo $INPUT| cut -d'_' -f 1)
mm=$(echo $INPUT| cut -d'_' -f 2)
yy=$(echo $INPUT| cut -d'_' -f 3 | tr -dc '0-9')
date=$(echo $dd'_'$mm'_'$yy)

### Loop through the dumps, extract .nc files, make sure to get any
### misplaced file where it belongs, update the .tar file, then "./diag"
for number in $(seq 0 $Dt_diag $a)
 do
 number_str=$( printf '%05d' $number )
  if [ -d "./t$number_str" ]
  then
  mv *3*$number_str.nc "t$number_str/"
  cd "t$number_str"

  rm -rf "t$number_str"
  tar -xvf *.tar
  mv *.nc t$number_str/
  rm -f *.tar
  tar -cf "Pack.tar" "t$number_str"
  mv t$number_str/*.nc .
  rm -rf t$number_str
  ../diag -t $number_str -d $date
  rm -f Mag3*
  rm -f Hsw3*
  rm -f Ele3*
  rm -f Vel3*
  rm -f Pre3*
  rm -f Den3*
  rm -f Atm3*

  rm -f Mom3*
#  rm -f p3*
  rm -f c3*
  rm -f b3*
  cd ..
  fi
 done

### Once the loop is over, place all the concatenated files together in the ncfiles/ folder
mv */Hsw*.nc .
mv */Magw*.nc .
mv */Thew*.nc .
mv */Elew*.nc .

mkdir ncfiles
mv *.nc ncfiles
rm -rf t0*

### Create images from .nc files
#cd ..
#cp -rf ~/bin/Python/ .
#cd Python/
#ipython plot_everything.py


