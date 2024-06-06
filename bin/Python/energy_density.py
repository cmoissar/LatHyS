# -*- coding: utf-8 -*-
import netCDF4 as nc
import pylab as pl
import scipy as scipy
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec
import gc
import sys

from class_data import *
from read_netcdf import *
from magfield_functions import *
from particle_functions import *
from plot_color_maps import *

def calc_dynamic_pressure(data_particle) :
    # Returns the dynamic pressure (in nPa)
    print('Calculate Pdyn')
    m_p = 1.67e-27                      # Mass of a proton
    data_particle.Pdyn = np.multiply(data_particle.n,
    np.multiply(data_particle.V,data_particle.V))*m_p*1.e21

    
def calc_mag_pressure(data_magfield,calc_grad=0) :
    # Returns the magnetic pressure (in nPa)
    print('Calculate Pmag')
    mu_0 = np.pi*4.e-7

    if calc_grad == 1:
        data_magfield.Pmag = np.multiply(data_magfield.B,data_magfield.B) \
        /(2.*mu_0)
    else :
        data_magfield.Pmag = np.multiply(data_magfield.B,data_magfield.B) \
        *1.e-9/(2.*mu_0)
        
        
def calc_thermal_pressure(data_particle,calc_grad=0) :
    # Returns the thermal pressure (in nPa)
    print('Calculate Pth')
    k_B = 1.3806488e-23
    # T is given in eV
    if calc_grad == 1 :
        data_particle.Pth = np.multiply(data_particle.n,data_particle.T)*k_B*1.e24*11605 
    else :
        data_particle.Pth = np.multiply(data_particle.n,data_particle.T)*k_B*1.e15*11605                               

def calc_total_pressure(data_magfield,data_particle) :
    # Returns the total pressure (in nPa)
    print('Calculate Ptot')
    data_magfield.Ptot = np.add(np.add(data_magfield.Pmag,
                                              data_particle.Pth),
                                data_particle.Pdyn)
    
def calc_beta(data_magfield,data_particle):
    # Returns the plasma beta, ratio between the plasma thermal pressure
    # and the magnetic pressure.

    data_magfield.beta = np.divide(data_particle.Pth,data_magfield.Pmag)

def calc_Mach_numbers(data_magfield,data_particle):
    # Calculates the Alfvenic, sonic and magnetosonic Mach numbers

    mu_0 = np.pi*4.e-7   # Permittivity of free space
    m_p  = 1.67e-27      # Mass of a proton (kg)
    k_B = 1.3806488e-23  # Boltzmann constant
    gamma = 5./3.        # Polytropic index

    # Alfven speed in km/s
    data_particle.Va = np.divide(data_magfield.B,
                                 np.sqrt(mu_0*m_p*data_particle.n))*1.e-15
    # Sound speed in km/s
    data_particle.cs = np.sqrt(gamma*k_B*data_particle.T*11605/m_p)*1.e-3

    # Magnetosonic speed in km/s
    data_particle.Vms = np.sqrt(np.add(np.multiply(data_particle.Va,data_particle.Va),
                                       np.multiply(data_particle.cs,data_particle.cs)))
    
    data_particle.Ma = np.divide(data_particle.V,data_particle.Va)
    data_particle.Ms = np.divide(data_particle.V,data_particle.cs)
    data_particle.Mms = np.divide(data_particle.V,data_particle.Vms)

if __name__ == '__main__' :
    
    pl.close()
    
    date = '20140101'
    date2 = '01_01_14'
    time = '00160'

    #filepath = '\\\\nas\\turc\\' + date + '\\'
    filepath = 'C:\\Users\\Lucile Turc\\Hybrid Simulations\\Outputs\\' + date + '\\'
    #filepath = 'D:\\Hybrid Simulations\\Outputs\\' + date + '\\'
    filepath_out = 'C:\\Users\\Lucile Turc\\Hybrid Simulations\\Plots\\' + date + '\\'
    #filepath = 'C:\\Users\\Lucile\\Work\\'#Nuages_magnétiques\\Python\\'
    #filepath_out = 'C:\\Users\\Lucile\\Work\\'#Nuages_magnétiques\\Python\\'
    file_Hsw = 'Hsw_'+date2+'_t' + time + '.nc'
    file_mag = 'Magw_'+date2+'_t' + time + '.nc'
    #file_ele = 'Elew_'+date2+'_t' + time + '.nc'

    XY = Plane(time)
    XZ = Plane(time)
    YZ = Plane(time)
    
    Hsw = Data(XY,XZ,YZ,time)
    Magw = Data(XY,XZ,YZ,time)
    #Elew = Data(XY,XZ,YZ,time)

    readNetcdfFile(filepath,file_Hsw,Hsw)
    readNetcdfFile(filepath,file_mag,Magw)
    #readNetcdfFile(filepath,file_ele,Elew,x_plane=-1)
    #plot_color_maps(Magw,'n',filepath_out,min_value=0.,max_value=15.)

    #sys.exit(0)
    Magw.x[:] = Magw.x[:]/Magw.c_omegapi
    Magw.y[:] = Magw.y[:]/Magw.c_omegapi
    Magw.z[:] = Magw.z[:]/Magw.c_omegapi

    Hsw.x[:] = Hsw.x[:]/Hsw.c_omegapi
    Hsw.y[:] = Hsw.y[:]/Hsw.c_omegapi
    Hsw.z[:] = Hsw.z[:]/Hsw.c_omegapi

##    Elew.x[:] = Elew.x[:]/Elew.c_omegapi
##    Elew.y[:] = Elew.y[:]/Elew.c_omegapi
##    Elew.z[:] = Elew.z[:]/Elew.c_omegapi
    
    calculate_velocity(Hsw.XY)
    calculate_velocity(Hsw.XZ)
    calculate_velocity(Hsw.YZ)

##    calculate_Etot(Elew.XY)
##    calculate_Etot(Elew.XZ)
##    calculate_Etot(Elew.YZ)

    #plot_color_maps(Elew,'Etot',filepath_out,min_value=0.,max_value=30.)
##    sys.exit(0)
    
    whereNaN = np.isnan(Hsw.XY.V)
    Hsw.XY.V[whereNaN] = 0.
    whereNaN = np.isnan(Hsw.XZ.V)
    Hsw.XZ.V[whereNaN] = 0.
    whereNaN = np.isnan(Hsw.YZ.V)
    Hsw.YZ.V[whereNaN] = 0.
    
##    plot_color_maps(Hsw,'V',filepath_out,min_value=0.,max_value=600.)
##
##    sys.exit(0)
    calculate_norm(Magw.XY)
    calculate_norm(Magw.XZ)
    calculate_norm(Magw.YZ)

    #print('n ',Hsw.XY.n[50,:])
    #print('B',Magw.XY.B.shape,Magw.XZ.B.shape,Magw.YZ.B.shape)

    calc_thermal_pressure(Hsw.XY)
    calc_thermal_pressure(Hsw.XZ)
    calc_thermal_pressure(Hsw.YZ)

    #sys.exit(0)
    
    calc_dynamic_pressure(Hsw.XY)
    calc_dynamic_pressure(Hsw.XZ)
    calc_dynamic_pressure(Hsw.YZ)

    whereNaN = np.isnan(Hsw.XY.Pdyn)
    Hsw.XY.Pdyn[whereNaN] = 0.
    whereNaN = np.isnan(Hsw.XZ.Pdyn)
    Hsw.XZ.Pdyn[whereNaN] = 0.
    whereNaN = np.isnan(Hsw.YZ.Pdyn)
    Hsw.YZ.Pdyn[whereNaN] = 0.
    
    #print('Pdyn ',Hsw.XY.Pdyn[50,:])
    calc_mag_pressure(Magw.XY)
    calc_mag_pressure(Magw.XZ)
    calc_mag_pressure(Magw.YZ)

    Msphere_XY = np.where(Hsw.XY.n < 1.5)
    Msphere_XZ = np.where(Hsw.XZ.n < 1.5)
    Msphere_YZ = np.where(Hsw.YZ.n < 1.5)

    Hsw.XY.Pth[Msphere_XY] = 0.
    Hsw.XZ.Pth[Msphere_XZ] = 0.
    Hsw.YZ.Pth[Msphere_YZ] = 0.

    Magw.XY.Pmag[Msphere_XY] = 0.
    Magw.XZ.Pmag[Msphere_XZ] = 0.
    Magw.YZ.Pmag[Msphere_YZ] = 0.

    calc_beta(Magw.XY,Hsw.XY)
    calc_beta(Magw.YZ,Hsw.YZ)
    calc_beta(Magw.XZ,Hsw.XZ)

    calc_Mach_numbers(Magw.XY,Hsw.XY)
    calc_Mach_numbers(Magw.YZ,Hsw.YZ)
    calc_Mach_numbers(Magw.XZ,Hsw.XZ)

    Hsw.XY.Ms[Msphere_XY] = 0.
    Hsw.XZ.Ms[Msphere_XZ] = 0.
    Hsw.YZ.Ms[Msphere_YZ] = 0.
    
    #print(Hsw.n_XY[10,:])

##    calc_electric_pressure(Elew.XY)
##    calc_electric_pressure(Elew.XZ)
##    calc_electric_pressure(Elew.YZ)
##        
##    calc_total_pressure(Magw.XY,Hsw.XY)
##    calc_total_pressure(Magw.XZ,Hsw.XZ)
##    calc_total_pressure(Magw.YZ,Hsw.YZ)
##    #print(Hsw.XY.Pdyn[50,:])
##    #print(Hsw.XY.Pdyn.shape)
##
##    print(Elew.XY.Ex[10,10],Elew.XY.Ey[10,10],Elew.XY.Ez[10,10])
##    print(Elew.XY.E[10,10])
##    
    #plot_color_maps(Magw,'B',filepath_out,min_value=0.,max_value=100.)
    #plot_color_maps(Magw,'Pmag',filepath_out,min_value=0,max_value=1,lang='En')
    #plot_color_maps(Hsw,'Pth',filepath_out,min_value=0,max_value=0.5,lang='En')
    #plot_color_maps(Hsw,'Pdyn',filepath_out,min_value=0,max_value=6,lang='En')
    #plot_color_maps(Elew,'Pel',filepath_out)#,min_value=0.,max_value=3.)
    #plot_color_maps(Magw,'Ptot',filepath_out,min_value=0.,max_value=5.,lang='En')
    plot_color_maps(Magw,'Beta',filepath_out,min_value=0.,max_value=5.,lang='En')
    #plot_color_maps(Magw,'Ma',filepath_out,min_value=0.,max_value=3.,lang='En')
    #plot_color_maps(Magw,'Ms',filepath_out,min_value=0.,max_value=12.,lang='En')
    #plot_color_maps(Hsw,'Mms',filepath_out,min_value=0.,max_value=3.,lang='En')
     
    del Hsw
    del Magw

