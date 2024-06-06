import netCDF4 as nc
import numpy as np
import os, re, glob, sys
from matplotlib import pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

####################    
### 3D functions ###
####################  

def import_3D(directory, tstep, verbose=False):
    tstep = f'{int(tstep):05d}'
    dictionnary_data = {}
    #Loop over the ncfiles to find 'Magw' and define useful params from it
    for nc_file in os.listdir(directory):
        if (not(tstep in nc_file) or not('Magw' in nc_file)):
    	    continue    

        f = nc.Dataset(directory + nc_file, 'r')
        print(f'Extracting simulation data for time {tstep}')
        if verbose:
            print(nc_file)
            print('Extracting general sinulation info from Magw...')
        
        coord_sys = f.variables['Coordinate_system'][:].data
        coord_name = b''
        for i in range(0, len(coord_sys)):
            coord_name = coord_name + coord_sys[i]

        coord_name = str(coord_name)[2:].split(' ')[0]

        if coord_name=='simu':
            sgn = -1
        else:
            sgn = +1

        dictionnary_data.update({'r_planet': np.float32(f.variables['r_planet'][:])})
        dictionnary_data.update({'c_omegapi': f.variables['phys_length'][:]})
        cwp = dictionnary_data['c_omegapi']
        dictionnary_data.update({'x': np.flip(f.variables['X_axis'][:])/cwp})
        dictionnary_data.update({'y': np.flip(f.variables['Y_axis'][:])/cwp})
        dictionnary_data.update({'z': f.variables['Z_axis'][:]/cwp})
        dictionnary_data.update({'gstep' : f.variables['gstep'][:]})

        global x,y,z
        x = dictionnary_data['x']
        y = dictionnary_data['y']
        z = dictionnary_data['z']
        nx,  ny,  nz  = len(x), len(y), len(z)
        nx0, ny0, nz0 = ( int(np.where(abs(x)==min(abs(x)))[0]),
                          int(np.where(abs(y)==min(abs(y)))[0]),
                          int(np.where(abs(z)==min(abs(z)))[0])  ) 

        #dictionnary_data.update({'x': x})
        #dictionnary_data.update({'y': y})
        #dictionnary_data.update({'z': z})
        dictionnary_data.update({'nx0': nx0})
        dictionnary_data.update({'ny0': ny0})
        dictionnary_data.update({'nz0': nz0})

        def extract_3D(A):
            '''
            This function extract the 3D data from the ncfile.
            Do not use this for large runs, where fooling around with the
            3D array is too much for most computers, and just takes ages for no valid reasons

            We use the array.transpose method for the following reason:
            arr = arr[z_axis,y_axis,x_axis]
            new_arr = arr.transpose(2, 1, 0)
            new_arr = new_arr[x_axis,y_axis,z_axis]
            '''  
            A = A.transpose(2, 1, 0)
            A = np.flip(A, axis=0)
            A = np.flip(A, axis=1)

            return A

    #Loop on the ncfiles a second time to extract 2D data (xy and xz planes)
    for nc_file in os.listdir(directory):
        if not(tstep in nc_file):
            continue    

        f = nc.Dataset(directory + nc_file, 'r')

        if 'Magw' in nc_file:
            if verbose:
                print('Reading Bx...')
            dictionnary_data.update({'Bx' : sgn*extract_3D(np.array(f.variables['Bx']))})
            if verbose:
                print('Reading By...')
            dictionnary_data.update({'By' : sgn*extract_3D(np.array(f.variables['By']))})
            if verbose:
                print('Reading Bz...')
            dictionnary_data.update({'Bz' : extract_3D(np.array(f.variables['Bz']))})

        elif 'Hsw' in nc_file:
            if verbose:
                print('Reading density...')
            dictionnary_data.update({'N' : extract_3D(np.array(f.variables['Density']))})
            if verbose:
                print('Reading Ux...')
            dictionnary_data.update({'Vx' : sgn*extract_3D(np.array(f.variables['Ux']))})
            if verbose:
                print('Reading Uy...')
            dictionnary_data.update({'Vy' : sgn*extract_3D(np.array(f.variables['Uy']))})               
            if verbose:
                print('Reading Uz...')
            dictionnary_data.update({'Vz' : extract_3D(np.array(f.variables['Uz']))})
            if verbose:
                print('Reading T...')
            dictionnary_data.update({'T' : extract_3D(np.array(f.variables['Temperature']))})

        elif 'Elew' in nc_file:
            if verbose:
                print('Reading Ex...')
            dictionnary_data.update({'Ex' : sgn*extract_3D(np.array(f.variables['Ex']))})
            if verbose:
                print('Reading Ey...')
            dictionnary_data.update({'Ey' : sgn*extract_3D(np.array(f.variables['Ey']))})
            if verbose:
                print('Reading Ez...')
            dictionnary_data.update({'Ez' : extract_3D(np.array(f.variables['Ez']))})

        f.close()
        
    return dictionnary_data

def load_data_3D(data_directory, t, verbose=False):
    simu_data = import_3D(data_directory, t, verbose)
    x, y, z = simu_data['x'], simu_data['y'], simu_data['z']
      
    N = simu_data['N']
    T = simu_data['T']
    Vx = simu_data['Vx']
    Vy = simu_data['Vy']
    Vz = simu_data['Vz']
    Ex = simu_data['Ex']
    Ey = simu_data['Ey']
    Ez = simu_data['Ez']
    Bx = simu_data['Bx']
    By = simu_data['By']
    Bz = simu_data['Bz']
    B = np.sqrt(simu_data['Bx']**2 + simu_data['By']**2 + simu_data['Bz']**2)
    E = np.sqrt(simu_data['Ex']**2 + simu_data['Ey']**2 + simu_data['Ez']**2)
    V = Vx**2+Vy**2+Vz**2
    nx0,  ny0,  nz0 = simu_data['nx0'], simu_data['ny0'], simu_data['nz0']
    
    return x, y, z, nx0, ny0, nz0, N, T, Vx, Vy, Vz, V, Ex, Ey, Ez, E, Bx, By, Bz, B


####################    
### 2D functions ###
####################  

def import_2D(directory, tstep, verbose=False):
    tstep = f'{int(tstep):05d}'
    dictionnary_data = {}
    #Loop over the ncfiles to find 'Magw' and define useful params from it
    for nc_file in os.listdir(directory):
        if (not(tstep in nc_file) or not('Magw' in nc_file)):
    	    continue    

        f = nc.Dataset(directory + nc_file, 'r')
        print(f'Extracting simulation data for time {tstep}')
        if verbose:
            print(nc_file)
            print('Extracting general sinulation info from Magw...')
        
        coord_sys = f.variables['Coordinate_system'][:].data
        coord_name = b''
        for i in range(0, len(coord_sys)):
            coord_name = coord_name + coord_sys[i]

        coord_name = str(coord_name)[2:].split(' ')[0]

        if coord_name=='simu':
            sgn = -1
        else:
            sgn = +1

        dictionnary_data.update({'r_planet': np.float32(f.variables['r_planet'][:])})
        dictionnary_data.update({'c_omegapi': f.variables['phys_length'][:]})
        cwp = dictionnary_data['c_omegapi']
        dictionnary_data.update({'x': np.flip(f.variables['X_axis'][:])/cwp})
        dictionnary_data.update({'y': np.flip(f.variables['Y_axis'][:])/cwp})
        dictionnary_data.update({'z': f.variables['Z_axis'][:]/cwp})
        dictionnary_data.update({'gstep' : f.variables['gstep'][:]})

        global x,y,z
        x = dictionnary_data['x']
        y = dictionnary_data['y']
        z = dictionnary_data['z']
        nx,  ny,  nz  = len(x), len(y), len(z)
        nx0, ny0, nz0 = ( int(np.where(abs(x)==min(abs(x)))[0]),
                          int(np.where(abs(y)==min(abs(y)))[0]),
                          int(np.where(abs(z)==min(abs(z)))[0])  ) 

        #Focus on a region of interest:
        slice_x = slice(-2*nx0-1,-1)   #this is for the un_flipped x
        slice_xy = (slice(nz0-1,nz0+1),slice(None),slice_x)
        slice_xz = (slice(None),slice(ny0-1,ny0+1),slice_x)

        x = x[:2*nx0]
        dictionnary_data.update({'x': x})
        dictionnary_data.update({'y': y})
        dictionnary_data.update({'z': z})
        dictionnary_data.update({'nx0': nx0})
        dictionnary_data.update({'ny0': ny0})
        dictionnary_data.update({'nz0': nz0})

        def extract_2D(A):
            '''
            This function extract the xy and xy plane data from the ncfile, without ever storing
            the whole 3D array. This is very important for big runs, where fooling around with the
            3D array is too much for most computers, and just takes ages for no valid reasons

            We use the array.transpose method for the following reason:
            arr = arr[z_axis,y_axis,x_axis]
            new_arr = arr.transpose(2, 1, 0)
            new_arr = new_arr[x_axis,y_axis,z_axis]
            '''  
            A_xy = np.average(A[slice_xy].transpose(2, 1, 0), 2)
            A_xy = np.flip(A_xy, axis=0)
            A_xy = np.flip(A_xy, axis=1)

            A_xz = np.average(A[slice_xz].transpose(2, 1, 0), 1)
            A_xz = np.flip(A_xz, axis=0)

            del(A)
            return A_xy, A_xz

    #Loop on the ncfiles a second time to extract 2D data (xy and xz planes)
    for nc_file in os.listdir(directory):
        if not(tstep in nc_file):
            continue    

        f = nc.Dataset(directory + nc_file, 'r')

        if 'Magw' in nc_file:
            if verbose:
                print('Reading Bx...')
            dictionnary_data.update({'Bx' : sgn*extract_2D(f.variables['Bx'])})
            if verbose:
                print('Reading By...')
            dictionnary_data.update({'By' : sgn*extract_2D(f.variables['By'])})
            if verbose:
                print('Reading Bz...')
            dictionnary_data.update({'Bz' : extract_2D(f.variables['Bz'])})

        elif 'Hsw' in nc_file:
            if verbose:
                print('Reading density...')
            dictionnary_data.update({'N' : extract_2D(f.variables['Density'])})
            if verbose:
                print('Reading Ux...')
            dictionnary_data.update({'Vx' : sgn*extract_2D(f.variables['Ux'])})
            if verbose:
                print('Reading Uy...')
            dictionnary_data.update({'Vy' : sgn*extract_2D(f.variables['Uy'])})               
            if verbose:
                print('Reading Uz...')
            dictionnary_data.update({'Vz' : extract_2D(f.variables['Uz'])})
            if verbose:
                print('Reading T...')
            dictionnary_data.update({'T' : extract_2D(f.variables['Temperature'])})

        elif 'Elew' in nc_file:
            if verbose:
                print('Reading Ex...')
            dictionnary_data.update({'Ex' : sgn*extract_2D(f.variables['Ex'])})
            if verbose:
                print('Reading Ey...')
            dictionnary_data.update({'Ey' : sgn*extract_2D(f.variables['Ey'])})
            if verbose:
                print('Reading Ez...')
            dictionnary_data.update({'Ez' : extract_2D(f.variables['Ez'])})

        f.close()
        
    return dictionnary_data

def load_data(data_directory, plane, t, verbose=False):
    simu_data = import_2D(data_directory, t, verbose)
    x, y, z = simu_data['x'], simu_data['y'], simu_data['z']
    if plane == 'xy':
        j = 0
        k = y
    if plane == 'xz':
        j = 1
        k = z
    
    Nxk = simu_data['N'][j]
    Vx_xk = simu_data['Vx'][j]
    Vy_xk = simu_data['Vy'][j]
    Vz_xk = simu_data['Vz'][j]
    Bxk = np.sqrt(simu_data['Bx'][j]**2 + simu_data['By'][j]**2 + simu_data['Bz'][j]**2)
    Exk = np.sqrt(simu_data['Ex'][j]**2 + simu_data['Ey'][j]**2 + simu_data['Ez'][j]**2)
    
    return x, k, Nxk, Vx_xk, Vy_xk, Vz_xk, Bxk, Exk

def plot_colormap(data_2D, plane):
    '''
    A small function to visualise things quickly
    '''
    if plane=='xy':
        p, k = 0, y
    if plane=='xz':
        p, k = 1, z

    A = data_2D

    plt.pcolor(x, k, A.T);ax=plt.gca(); ax.set_aspect('equal'); ax.invert_xaxis(); plt.show()
    
def colormap(x, k, Nxk, plane, tstep, k_shock_list=None, x_shock_list=None, k_cut=None, xlims=None, klims=None):

    # COLORMAP #
    plt.rcParams["figure.figsize"] = (8,8)    
    plt.pcolor(x, k, Nxk.T)
    if xlims==None:
        plt.xlim((max(x), min(x)))
    else:
        plt.xlim(xlims)
    if klims==None and not(k_shock_list==None):
        plt.ylim((min(k_shock_list), max(k_shock_list)))
    else:
        plt.ylim(klims)
    plt.ylabel(f'{plane[1]} (di)', fontsize=18)
    plt.xlabel('x (di)', fontsize=18)
    ax=plt.gca()
    ax.set_aspect('equal');

    # TIME STEP #
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    # print(f"xmax = {xmax}, xmin = {xmin}, ymax = {ymax}, ymin = {ymin}")
    props = dict(boxstyle='round', facecolor='wheat', alpha=1)
    plt.text(xmin+(xmax-xmin)*0.05, ymax-(ymax-ymin)*0.05, f'time {tstep}', fontsize=12, bbox=props)

    # SHOCK #
    if not(x_shock_list==None) and not(k_shock_list==None):
        ax.scatter(x_shock_list, k_shock_list, color='red',s=2.5)

    # hline #
    if k_cut:
        ax.hlines(k_cut, xmin, xmax, color='red', linewidth=2.5)
        
    # COLORBAR #
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    plt.clim(0,40)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(cax=cax)
    cbar.set_label(label=r'N(cm$^{-3}$)', rotation=270, fontsize = 16, weight="bold", labelpad=20)

    Output_directory = './analysis/'
    plt.savefig(Output_directory + 'shock_detection' + plane + tstep + '.png', bbox_inches='tight', pad_inches=0.05)
    plt.show()
    
####################    
### 1D functions ###
####################  

def import_1D(directory, tstep, y_cut=0, z_cut=0, verbose=False):
    tstep = f'{int(tstep):05d}'
    dictionnary_data = {}
    #Loop over the ncfiles to find 'Magw' and define useful params from it
    for nc_file in os.listdir(directory):
        if (not(tstep in nc_file) or not('Magw' in nc_file)):
            continue    

        f = nc.Dataset(directory + nc_file, 'r')
        print(f'Extracting simulation data for time {tstep}')
        if verbose:
            print(nc_file)
            print('Extracting general sinulation info from Magw...')
        
        coord_sys = f.variables['Coordinate_system'][:].data
        coord_name = b''
        for i in range(0, len(coord_sys)):
            coord_name = coord_name + coord_sys[i]

        coord_name = str(coord_name)[2:].split(' ')[0]

        if coord_name=='simu':
            sgn = -1
        else:
            sgn = +1

        dictionnary_data.update({'r_planet': np.float32(f.variables['r_planet'][:])})
        dictionnary_data.update({'c_omegapi': f.variables['phys_length'][:]})
        cwp = dictionnary_data['c_omegapi']
        dictionnary_data.update({'x': np.flip(f.variables['X_axis'][:])/cwp})
        dictionnary_data.update({'y': np.flip(f.variables['Y_axis'][:])/cwp})
        dictionnary_data.update({'z': f.variables['Z_axis'][:]/cwp})
        dictionnary_data.update({'gstep' : f.variables['gstep'][:]})

        global x,y,z
        x = dictionnary_data['x']
        y = dictionnary_data['y']
        z = dictionnary_data['z']
        nx,  ny,  nz  = len(x), len(y), len(z)
        nx0, ny0, nz0 = ( int(np.where(abs(x)==min(abs(x)))[0]),
                          int(np.where(abs(y-y_cut)==min(abs(y-y_cut)))[0]),
                          int(np.where(abs(z-z_cut)==min(abs(z-z_cut)))[0])  ) 

        dictionnary_data.update({'x': x})
        dictionnary_data.update({'y': y})
        dictionnary_data.update({'z': z})
        dictionnary_data.update({'nx0': nx0})
        dictionnary_data.update({'ny0': ny0})
        dictionnary_data.update({'nz0': nz0})

        def extract_1D(A):
            '''
            This function extract the line data at y=z=0 from the ncfile, without ever storing
            the whole 3D array. This is very important for big runs, where fooling around with the
            3D array is too much for most computers, and just takes ages for no valid reasons.

            We use the array.transpose method for the following reason:
            arr = arr[z_axis,y_axis,x_axis]
            new_arr = arr.transpose(2, 1, 0)
            new_arr = new_arr[x_axis,y_axis,z_axis]
            '''  
            A_line = A[nz0, ny0, :]
            A_line = np.flip(A_line, axis=0)

            del(A)
            return A_line

    #Loop on the ncfiles a second time to extract 1D data 
    for nc_file in os.listdir(directory):
        if not(tstep in nc_file):
            continue    

        f = nc.Dataset(directory + nc_file, 'r')

        if 'Magw' in nc_file:
            if verbose:
                print('Reading Bx...')
            dictionnary_data.update({'Bx' : sgn*extract_1D(f.variables['Bx'])})
            if verbose:
                print('Reading By...')
            dictionnary_data.update({'By' : sgn*extract_1D(f.variables['By'])})
            if verbose:
                print('Reading Bz...')
            dictionnary_data.update({'Bz' : extract_1D(f.variables['Bz'])})

        elif 'Hsw' in nc_file:
            if verbose:
                print('Reading density...')
            dictionnary_data.update({'N' : extract_1D(f.variables['Density'])})
            if verbose:
                print('Reading Ux...')
            dictionnary_data.update({'Vx' : sgn*extract_1D(f.variables['Ux'])})
            if verbose:
                print('Reading Uy...')
            dictionnary_data.update({'Vy' : sgn*extract_1D(f.variables['Uy'])})               
            if verbose:
                print('Reading Uz...')
            dictionnary_data.update({'Vz' : extract_1D(f.variables['Uz'])})
            if verbose:
                print('Reading T...')
            dictionnary_data.update({'T' : extract_1D(f.variables['Temperature'])})

        elif 'Elew' in nc_file:
            if verbose:
                print('Reading Ex...')
            dictionnary_data.update({'Ex' : sgn*extract_1D(f.variables['Ex'])})
            if verbose:
                print('Reading Ey...')
            dictionnary_data.update({'Ey' : sgn*extract_1D(f.variables['Ey'])})
            if verbose:
                print('Reading Ez...')
            dictionnary_data.update({'Ez' : extract_1D(f.variables['Ez'])})

        f.close()
        
    return dictionnary_data

def load_data_1D(data_directory, t, y_cut=0, z_cut=0, verbose=False):
    simu_data = import_1D(data_directory, t, y_cut=y_cut, z_cut=z_cut, verbose=verbose)
    x, y, z = simu_data['x'], simu_data['y'], simu_data['z']
    
    N = simu_data['N']
    Vx = simu_data['Vx']
    Vy = simu_data['Vy']
    Vz = simu_data['Vz']
    B = np.sqrt(simu_data['Bx']**2 + simu_data['By']**2 + simu_data['Bz']**2)
    E = np.sqrt(simu_data['Ex']**2 + simu_data['Ey']**2 + simu_data['Ez']**2)
    
    return x, N, Vx, Vy, Vz, B, E
