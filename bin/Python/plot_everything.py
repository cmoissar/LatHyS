# -*- coding: utf-8 -*-
import pylab as pl
import numpy as np

from class_data import Plane, Data
from read_netcdf import readNetcdfFile
from particle_functions import calculate_velocity
from plot_color_maps import plot_color_maps
import os
import glob
import re

from energy_density import *       #contains all the calc_xxxxx functions

#import pdb
     #   pdb.set_trace()



if __name__ == '__main__' :
     

        ##############################################################
        # Modify filepath and filepath_out (depends on the computer) #
        ##############################################################
    
    Cluster = 'Zoidberg'
    Cluster = 'Curie'
    Cluster = 'Occ'
    
    os.system('ulimit -s 32768')
 
    liste = ['.']
 
    for run_name in liste :

        filepath = '../ncfiles/'

        filepath_out = filepath + '..' + '/Images/'
        
        if not os.path.exists(filepath_out):
            os.mkdir(filepath_out)
            
            
            ##############################################################
            #         Retrieve useful parameters from the launcher       #
            ##############################################################
                    

        param_list = ['TMAX', 'DT']
        file = filepath + '../' + 'Code/Soumission/' + 'launcher.sh'
        
        # read the file into a list of lines
        with open(file,'r', encoding="ISO-8859-1") as f:
 
            lines = f.read().split("\n")

        for param in param_list:

            # iterate over lines, and print out line numbers which contain
            # the word of interest.
            for i, line in enumerate(lines):
                if "echo" in line:
                    continue
                if param in line:  # or word in line.split() to search for full words
                    n_line = i
                    break            
            exec(open(file,encoding = "ISO-8859-1").readlines()[n_line])        
               
        #NHM,DT and NTREG are implicity defined from the previous "for loop"    
        tmax  = TMAX
        NHM = tmax / DT
        tstep = 1

            # In case something went wrong, manually enter the values
        # tmax  = 250
        # tstep = 10

        # Date of the simulation, as written in
        # the ncfiles.
   
        date = re.search('w_(.+?)_t', glob.glob(filepath+'Hsw*_t00000*')[0]).group(1)
        
        for time in range(0,tmax+1,tstep):

            try :

                time = '%05d' % time    # Time of the output

                pl.close()

                file_Hsw = 'Hsw_'+date+'_t' + time + '.nc'
                file_mag = 'Magw_'+date+'_t' + time + '.nc'

                XY = Plane(time)
                XZ = Plane(time)
                YZ = Plane(time)

                Hsw = Data(XY,XZ,YZ,time)
                Magw = Data(XY,XZ,YZ,time)

                # Read the data in the netcdf file in 3 planes at x = x_plane = constant,
                # y = y_plane = constant and z = z_plane = constant
                # Default: x_plane = 0 : terminator plane
                #          y_plane = 0 : noon-midnight meridian plane
                #          z_plane = 0 : equatorial plane
                # (if the obstacle is at the centre of the simulation domain)
                # Also: x_plane = -1 : exit face of the simulation box
                # Also: x_plane = 2/3 : Planet near the exit face but not stuck to it
                print('Reading file...')
                readNetcdfFile(filepath,file_Hsw,Hsw, x_plane=0)
                readNetcdfFile(filepath,file_mag,Magw, x_plane=0)

                # Axes in normalised units
                Hsw.x[:] = Hsw.x[:]/Hsw.c_omegapi
                Hsw.y[:] = Hsw.y[:]/Hsw.c_omegapi
                Hsw.z[:] = Hsw.z[:]/Hsw.c_omegapi
                # Axes in normalised units
                Magw.x[:] = Magw.x[:]/Magw.c_omegapi
                Magw.y[:] = Magw.y[:]/Magw.c_omegapi
                Magw.z[:] = Magw.z[:]/Magw.c_omegapi

                # Calculate the total velocity
                calculate_velocity(Hsw.XY)
                calculate_velocity(Hsw.XZ)
                calculate_velocity(Hsw.YZ)
                # Calculate the thermal pressure
                calc_thermal_pressure(Hsw.XY)
                calc_thermal_pressure(Hsw.XZ)
                calc_thermal_pressure(Hsw.YZ)
                # Calculate the magnetic field magnitude
                calculate_norm(Magw.XY)
                calculate_norm(Magw.XZ)
                calculate_norm(Magw.YZ)
                # Calculate the magnetic pressure
                calc_mag_pressure(Magw.XY)
                calc_mag_pressure(Magw.XZ)
                calc_mag_pressure(Magw.YZ)

                # We define the magnetosphere as the region where the density is below 1.5/cc
                Msphere_XY = np.where(Hsw.XY.n < 1.5)
                Msphere_XZ = np.where(Hsw.XZ.n < 1.5)
                Msphere_YZ = np.where(Hsw.YZ.n < 1.5)

                # We put to 0 the thermal pressure and the magnetic pressure
                # everywhere inside the magnetosphere as defined above
                Hsw.XY.Pth[Msphere_XY] = 0.
                Hsw.XZ.Pth[Msphere_XZ] = 0.
                Hsw.YZ.Pth[Msphere_YZ] = 0.

                Magw.XY.Pmag[Msphere_XY] = 0.
                Magw.XZ.Pmag[Msphere_XZ] = 0.
                Magw.YZ.Pmag[Msphere_YZ] = 0.

                # Calculation of the plasma beta
                calc_beta(Magw.XY,Hsw.XY)
                calc_beta(Magw.YZ,Hsw.YZ)
                calc_beta(Magw.XZ,Hsw.XZ)

                # Plot the velocity in the 3 planes
                # min_value and max_value are the lower and upper limits of the color scale
                # if they are both set to 0, color scale is automatic

                print('Plot velocity color map...')
                plot_color_maps(Hsw,'V',filepath_out+date,min_value=0,max_value=0)

                print('Plot density color map...')
                plot_color_maps(Hsw,'n',filepath_out+date,min_value=0,max_value=0)

                print('Plot magnetic field color map...')
                plot_color_maps(Magw,'B',filepath_out+date,min_value=0,max_value=0)

                # Plot the plasma beta in the 3 planes
                # min_value and max_value are the lower and upper limits of the color scale
                plot_color_maps(Magw,'Beta',filepath_out+date,min_value=0.,max_value=1.,lang='En')

                del Hsw
                del Magw

            except :
                pass
