# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')    #choose backend before anything else, so that 
                         #the code doesn't try to open any window 
                         #(cannot do it on cluster)

import pylab as pl
import numpy as np

from class_data import Plane, Data
from read_netcdf import readNetcdfFile
from plot_color_maps import plot_color_maps
import os
import glob
import re
from pathlib2 import Path

from energy_density import *       #contains all the calc_xxxxx functions

#import pdb
     #   pdb.set_trace()



if __name__ == '__main__' :

        ##############################################################
        # Modify filepath and filepath_out (depends on the computer) #
        ##############################################################
    
    liste = [ '../ncfiles/' ]

    for filepath in liste : 

        filepath_out = filepath + '..' + '/Images/'
        
        if not os.path.exists(filepath_out):
            os.mkdir(filepath_out)
            
            
            ##############################################################
            #         Retrieve useful parameters from the launcher       #
            ##############################################################
                    
        param_list = ['TMAX', 'DT']
        file = filepath + '../Code/Soumission/launcher.sh'
        
        # read the file into a list of lines
        with open(file) as f:
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
            exec(open(file).readlines()[n_line])        
               
        #NHM,DT and NTREG are implicity defined from the previous "for loop"    
        tmax  = TMAX
        NHM   = int(tmax / DT)
        tstep = 1

        # Retrieve Date of the simulation, as written in
        # the ncfiles
        date = re.search('w_(.+?)_t', glob.glob(filepath + 'Hsw*_t00000*')[0]).group(1)
             
        # Set-up the x_axis (otherwise it changes overtime because Hsw is weird
        time0 = '00000'
        file_mag = 'Magw_'+date+'_t' + time0 + '.nc'
        XY = Plane(time0)
        XZ = Plane(time0)
        YZ = Plane(time0)
        Magw0 = Data(XY,XZ,YZ,time0)
        readNetcdfFile(filepath,file_mag,Magw0,x_plane=-1)
        # Axes in normalised units
        x0 = Magw0.x[:]/Magw0.c_omegapi

        #for time in List :
        for time in range (0,tmax+1,tstep):
            
            time = '%05d' % time    # Time of the output
        
            pl.close()    
    
            file_Hsw = 'Hsw_' +date+'_t' + time + '.nc'
            file_mag = 'Magw_'+date+'_t' + time + '.nc'
            file_ele = 'Elew_'+date+'_t' + time + '.nc'

            f_Hsw = Path(filepath + file_Hsw)
            f_mag = Path(filepath + file_mag)

            if not (f_Hsw.is_file() and f_mag.is_file()):
                print("files are missing for timestep: " + str(time))
                continue

    
            XY = Plane(time)
            XZ = Plane(time)
            YZ = Plane(time)
        
            Hsw  = Data(XY,XZ,YZ,time)
            Magw = Data(XY,XZ,YZ,time)
            Elew = Data(XY,XZ,YZ,time)
            

    
            # Read the data in the netcdf file in 3 planes at x = x_plane = constant,
            # y = y_plane = constant and z = z_plane = constant
            # Default: x_plane = 0 : terminator plane
            #          y_plane = 0 : noon-midnight meridian plane
            #          z_plane = 0 : equatorial plane
            # (if the obstacle is at the centre of the simulation domain)
            # Also: x_plane = -1 : exit face of the simulation box
            print('Reading file...') 
            readNetcdfFile(filepath,file_Hsw,Hsw,x_plane=-1)
            readNetcdfFile(filepath,file_mag,Magw,x_plane=-1)
            readNetcdfFile(filepath,file_ele,Elew,x_plane=-1)

            # import pdb
            # pdb.set_trace()
            
            # Axes in normalised units
            Hsw.x[:] = Hsw.x[:]/Hsw.c_omegapi
            Hsw.y[:] = Hsw.y[:]/Hsw.c_omegapi
            Hsw.z[:] = Hsw.z[:]/Hsw.c_omegapi
            # Axes in normalised units
            Magw.x[:] = Magw.x[:]/Magw.c_omegapi
            Magw.y[:] = Magw.y[:]/Magw.c_omegapi
            Magw.z[:] = Magw.z[:]/Magw.c_omegapi
            # Axes in normalised units
            Elew.x[:] = Elew.x[:]/Elew.c_omegapi
            Elew.y[:] = Elew.y[:]/Elew.c_omegapi
            Elew.z[:] = Elew.z[:]/Elew.c_omegapi
            
            # Calculate the total velocity
            from particle_functions import calculate_velocity
            calculate_velocity(Hsw.XY)
            calculate_velocity(Hsw.XZ)
            calculate_velocity(Hsw.YZ)
            # Calculate the thermal pressure
            from particle_functions import calc_thermal_pressure
            calc_thermal_pressure(Hsw.XY)
            calc_thermal_pressure(Hsw.XZ)
            calc_thermal_pressure(Hsw.YZ)
            # Calculate the electric field magnitude
            from elefield_functions import calculate_norm
            calculate_norm(Elew.XY)
            calculate_norm(Elew.XZ)
            calculate_norm(Elew.YZ)
            # Calculate the magnetic field magnitude
            from magfield_functions import calculate_norm
            calculate_norm(Magw.XY)
            calculate_norm(Magw.XZ)
            calculate_norm(Magw.YZ)
            # Calculate the magnetic pressure
            from energy_density import calc_mag_pressure
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

            import matplotlib.gridspec as gridspec
            def cm2inch(value):
                return value / 2.54
            plt.figure(figsize=(cm2inch(21), cm2inch(29.7)))
            nb_lines = 8
            gs = gridspec.GridSpec(nb_lines, 1)
            fsize = 12

            def line_plot(name,unit,cl_data,data_XY,data_XZ,i, label=None):
                axe = plt.subplot(gs[i])
                # moyenne arithmétique un peu bête. Ne prend pas en compte la taille en Y ou Z
                signal = (np.average(data_XY[:-1, :], 0) + np.average(data_XZ[:-1, :], 0)) / 2
                plt.plot(x0, signal, label=label)
                axe.set_xlim([np.min(x0), np.max(x0)])
                plt.gca().invert_xaxis()
                axe.get_xaxis().set_visible(False)
                plt.title(name, fontsize=fsize)
                plt.ylabel(unit, fontsize=fsize)
                plt.legend(loc='best')
                if (i == nb_lines-1):
                    axe.get_xaxis().set_visible(True)
                    plt.xlabel('x', fontsize=fsize)


            flatui = ["#0082c8"	, "#3cb44b", "#808080", "#e74c3c", "#34495e", "#2ecc71"]

            plt.rcParams['axes.prop_cycle'] = plt.cycler(color=flatui)

            line_plot("magnetic field"       , "$(nT)$", Magw, Magw.XY.B , Magw.XZ.B , 0)
            line_plot("magnetic field components", "$(nT)$", Magw, Magw.XY.Bx, Magw.XZ.Bx, 1, label="Bx")
            line_plot("magnetic field components", "$(nT)$", Magw, Magw.XY.By, Magw.XZ.By, 1, label="By")
            line_plot("magnetic field components", "$(nT)$", Magw, Magw.XY.Bz, Magw.XZ.Bz, 1, label="Bz")

            line_plot("electric field components", "$(V/m)$", Elew, Elew.XY.Ex, Elew.XZ.Ex, 2, label='Ex')
            line_plot("electric field components", "$(V/m)$", Elew, Elew.XY.Ey, Elew.XZ.Ey, 2, label='Ey')
            line_plot("electric field components", "$(V/m)$", Elew, Elew.XY.Ez, Elew.XZ.Ez, 2, label='Ez')

            line_plot("ion velocity"    , "$(km/s)$"    , Hsw, Hsw.XY.V, Hsw.XZ.V, 3)
            line_plot("ion velocity components"    , "$(km/s)$"    , Hsw, Hsw.XY.Vx, Hsw.XZ.Vx, 4, label='Vx')
            line_plot("ion velocity components"    , "$(km/s)$"    , Hsw, Hsw.XY.Vy, Hsw.XZ.Vy, 4, label='Vy')
            line_plot("ion velocity components"    , "$(km/s)$"    , Hsw, Hsw.XY.Vz, Hsw.XZ.Vz, 4, label='Vz')
            line_plot("ion density"     , "$(nb/cm^3)$" , Hsw, Hsw.XY.n, Hsw.XZ.n, 5)
            line_plot("ion temperature" , "$(eV)$"      , Hsw, Hsw.XY.T, Hsw.XZ.T, 6)

            line_plot("Beta", " ",  Magw, Magw.XY.beta, Magw.XZ.beta, 7)

            plt.suptitle(' t = '+ Hsw.time +
              r"\boldmath{$\, \Omega_{\mathrm{ci}}^{-1}$}",
                 fontsize=fsize*1.5,y=0.92,weight="bold" )

                        # Fine-tune Page; make subplots farther from each other.
            plt.subplots_adjust(hspace=cm2inch(0.7))
            plt.subplots_adjust(wspace=0)

            fileout = date + 'lines' + 't=' + Hsw.time + '.png'

            print(filepath_out + fileout)
            plt.savefig(filepath_out + fileout, dpi=150)
            print('Figure saved as ', filepath_out + fileout)


            del Hsw
            del Magw
            del Elew
