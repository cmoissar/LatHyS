import netCDF4 as nc
from class_data import *
import sys
import numpy as np

def readNetcdfFile(filepath,filename,str_data,x_plane=0,y_plane=0,z_plane=0, plane_sup=0,plane_inf=0) :

    f = nc.Dataset(filepath+filename,'r')

    str_data.c_omegapi = f.variables['phys_length'][:]
    str_data.x  = f.variables['X_axis'][:]
    str_data.y  = f.variables['Y_axis'][:]
    str_data.z  = f.variables['Z_axis'][:]

    str_data.gstep = f.variables['gstep'][:]

    coord_sys = f.variables['Coordinate_system'][:].data

    coord_name = b''
    for i in range(0, len(coord_sys)):
        coord_name = coord_name + coord_sys[i]
    coord_name = str(coord_name)[2:].split(' ')[0]

    if coord_name=='simu':
        sgn = -1
    else:
        sgn = +1
 
    if not(x_plane == 0) :
        x_plane = int(str_data.x.shape[0]*x_plane)
    if x_plane == 0 :
        x_plane = np.where(abs(str_data.x)<str_data.c_omegapi*str_data.gstep[0])[0][0] # int(str_data.x.shape[0]/2)
    if x_plane == -1 :
        x_plane = str_data.x.shape[0] - 5
    if y_plane == 0 :
        y_plane = int(str_data.y.shape[0]/2)
    if z_plane == 0 :
        z_plane = int(str_data.y.shape[0]/2)
    if plane_sup == 1:
        x_plane = x_plane - 1
        y_plane = y_plane - 1
        z_plane = z_plane + 1
    elif plane_inf == 1:
        x_plane = x_plane + 1
        y_plane = y_plane + 1
        z_plane = z_plane - 1

    if filename[0:4] == 'Magw' :
        print('Reading Bx...')
        str_data.XY.Bx = f.variables['Bx'][z_plane,:,:]
        str_data.XZ.Bx = f.variables['Bx'][:,y_plane,:]
        str_data.YZ.Bx = f.variables['Bx'][:,:,x_plane]


        print('Reading By...')
        str_data.XY.By = f.variables['By'][z_plane,:,:]
        str_data.XZ.By = f.variables['By'][:,y_plane,:]
        str_data.YZ.By = f.variables['By'][:,:,x_plane]

        print('Reading Bz...')
        str_data.XY.Bz = f.variables['Bz'][z_plane,:,:]
        str_data.XZ.Bz = f.variables['Bz'][:,y_plane,:]
        str_data.YZ.Bz = f.variables['Bz'][:,:,x_plane]      
   
    elif filename[0:3] == 'Hsw' :            
        #print(filename) 
        #print(str_data.n_XY.shape)
        #print(f.variables['Density'].shape)
        #print(f.variables['Density'][10,10,10])
        print('Reading density...')
        str_data.XY.n = f.variables['Density'][z_plane,:,:]
        str_data.XZ.n = f.variables['Density'][:,y_plane,:]
        str_data.YZ.n = f.variables['Density'][:,:,x_plane]

        print('Reading Ux...')
        str_data.XY.Vx = f.variables['Ux'][z_plane,:,:]
        str_data.XZ.Vx = f.variables['Ux'][:,y_plane,:]
        str_data.YZ.Vx = f.variables['Ux'][:,:,x_plane]

        print('Reading Uy...')
        str_data.XY.Vy = f.variables['Uy'][z_plane,:,:]
        str_data.XZ.Vy = f.variables['Uy'][:,y_plane,:]
        str_data.YZ.Vy = f.variables['Uy'][:,:,x_plane]

        print('Reading Uz...')                  
        str_data.XY.Vz = f.variables['Uz'][z_plane,:,:]
        str_data.XZ.Vz = f.variables['Uz'][:,y_plane,:]
        str_data.YZ.Vz = f.variables['Uz'][:,:,x_plane]                  

        print('Reading T...')
        str_data.XY.T = f.variables['Temperature'][z_plane,:,:]
        str_data.XZ.T = f.variables['Temperature'][:,y_plane,:]
        str_data.YZ.T = f.variables['Temperature'][:,:,x_plane]

    elif filename[0:4] == 'Elew':

        print('Reading Ex...')
        str_data.XY.Ex = f.variables['Ex'][z_plane,:,:]*1e+6
        str_data.XZ.Ex = f.variables['Ex'][:,y_plane,:]*1e+6
        str_data.YZ.Ex = f.variables['Ex'][:,:,x_plane]*1e+6

        print('Reading Ey...')
        str_data.XY.Ey = f.variables['Ey'][z_plane,:,:]*1e+6
        str_data.XZ.Ey = f.variables['Ey'][:,y_plane,:]*1e+6
        str_data.YZ.Ey = f.variables['Ey'][:,:,x_plane]*1e+6

        print('Reading Ez...')
        str_data.XY.Ez = f.variables['Ez'][z_plane,:,:]*1e+6
        str_data.XZ.Ez = f.variables['Ez'][:,y_plane,:]*1e+6
        str_data.YZ.Ez = f.variables['Ez'][:,:,x_plane]*1e+6
        
    print('Close file and return...')              
    f.close()
    return
    
def readNetcdfFile_point(filepath,filename,str_data,coord,ifile) :

    f = nc.Dataset(filepath+filename,'r')
    #for v in f.variables:
    #    print(v)

    if filename[0:4] == 'Magw' :

        str_data.Bx[ifile] = f.variables['Bx'][coord[2],coord[1],coord[0]] 
        str_data.By[ifile] = f.variables['By'][coord[2],coord[1],coord[0]] 
        str_data.Bz[ifile] = f.variables['Bz'][coord[2],coord[1],coord[0]]  

    elif filename[0:3] == 'Hsw' :            

        str_data.n[ifile] = f.variables['Density'][coord[2],coord[1],coord[0]]
        str_data.Vx[ifile] = f.variables['Ux'][coord[2],coord[1],coord[0]]
        str_data.Vy[ifile] = f.variables['Uy'][coord[2],coord[1],coord[0]]
        str_data.Vz[ifile] = f.variables['Uz'][coord[2],coord[1],coord[0]]
        str_data.T[ifile] = f.variables['Temperature'][coord[2],coord[1],coord[0]]              
      
    f.close()
    return

def readNetcdfFile_point_average(filepath,filename,str_data,coord,ifile) :

    f = nc.Dataset(filepath+filename,'r')
    #for v in f.variables:
    #    print(v)

    if filename[0:4] == 'Magw' :

        str_data.Bx[ifile] = (f.variables['Bx'][coord[2],coord[1],coord[0]] + \
                             f.variables['Bx'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['Bx'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['Bx'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['Bx'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['Bx'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['Bx'][coord[2],coord[1],coord[0]-1])/7.
        
        str_data.By[ifile] = (f.variables['By'][coord[2],coord[1],coord[0]] + \
                             f.variables['By'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['By'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['By'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['By'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['By'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['By'][coord[2],coord[1],coord[0]-1])/7.
     
        str_data.Bz[ifile] = (f.variables['Bz'][coord[2],coord[1],coord[0]] + \
                             f.variables['Bz'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['Bz'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['Bz'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['Bz'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['Bz'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['Bz'][coord[2],coord[1],coord[0]-1])/7.
     
    elif filename[0:3] == 'Hsw' :
        
        str_data.n[ifile] = (f.variables['Density'][coord[2],coord[1],coord[0]] + \
                             f.variables['Density'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['Density'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['Density'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['Density'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['Density'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['Density'][coord[2],coord[1],coord[0]-1])/7.
    
        str_data.Vx[ifile] = (f.variables['Ux'][coord[2],coord[1],coord[0]] + \
                             f.variables['Ux'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['Ux'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['Ux'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['Ux'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['Ux'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['Ux'][coord[2],coord[1],coord[0]-1])/7.
    
        str_data.Vy[ifile] = (f.variables['Uy'][coord[2],coord[1],coord[0]] + \
                             f.variables['Uy'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['Uy'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['Uy'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['Ux'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['Uy'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['Uy'][coord[2],coord[1],coord[0]-1])/7.

        str_data.Vz[ifile] = (f.variables['Uz'][coord[2],coord[1],coord[0]] + \
                             f.variables['Uz'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['Uz'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['Uz'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['Uz'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['Uz'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['Uz'][coord[2],coord[1],coord[0]-1])/7.

        str_data.T[ifile] = (f.variables['Temperature'][coord[2],coord[1],coord[0]] + \
                             f.variables['Temperature'][coord[2]+1,coord[1],coord[0]] + \
                             f.variables['Temperature'][coord[2]-1,coord[1],coord[0]] + \
                             f.variables['Temperature'][coord[2],coord[1]+1,coord[0]] + \
                             f.variables['Temperature'][coord[2],coord[1]-1,coord[0]] + \
                             f.variables['Temperature'][coord[2],coord[1],coord[0]+1] + \
                             f.variables['Temperature'][coord[2],coord[1],coord[0]-1])/7.

    f.close()
    return


def readNetcdfFile_line(filepath,filename,str_data,xline=0,yline=0,zline=0) :

    f = nc.Dataset(filepath+filename,'r')

##    for v in f.variables:
##        print(v)
 
    
    
    str_data.ref_dens = f.variables['phys_density'][:]
    str_data.c_omegapi = f.variables['phys_length'][:]
    #print(xline,yline,zline)
    if xline==-1:

        str_data.x = f.variables['X_axis'][:]
        
        if filename[0:4] == 'Magw' :

            str_data.Bx = f.variables['Bx'][zline,yline,:] 
            str_data.By = f.variables['By'][zline,yline,:]  
            str_data.Bz = f.variables['Bz'][zline,yline,:]  
        elif filename[0:3] == 'Hsw':
            
            str_data.n = f.variables['Density'][zline,yline,:] 
            str_data.Vx = f.variables['Ux'][zline,yline,:]
            str_data.Vy = f.variables['Uy'][zline,yline,:]
            str_data.Vz = f.variables['Uz'][zline,yline,:]
            str_data.T = f.variables['Temperature'][zline,yline,:]

    elif yline == -1:

        str_data.y = f.variables['Y_axis'][:]
        
        if filename[0:4] == 'Magw' :

            str_data.Bx = f.variables['Bx'][zline,:,xline] 
            str_data.By = f.variables['By'][zline,:,xline] 
            str_data.Bz = f.variables['Bz'][zline,:,xline] 
        elif filename[0:3] == 'Hsw':
            #print(xline,zline)
            #print(f.variables['Density'].shape)
            #print(f.variables['Density'][zline,yline,100])
            str_data.n = f.variables['Density'][zline,:,xline] 
            str_data.Vx = f.variables['Ux'][zline,:,xline] 
            str_data.Vy = f.variables['Uy'][zline,:,xline] 
            str_data.Vz = f.variables['Uz'][zline,:,xline] 
            str_data.T = f.variables['Temperature'][zline,:,xline] 

    elif zline == -1:

        str_data.z = f.variables['Z_axis'][:]
        
        if filename[0:4] == 'Magw' :

            str_data.Bx = f.variables['Bx'][:,yline,xline] 
            str_data.By = f.variables['By'][:,yline,xline]
            str_data.Bz = f.variables['Bz'][:,yline,xline]
        elif filename[0:3] == 'Hsw':
            #print(xline,zline)
            #print(f.variables['Density'].shape)
            #print(f.variables['Density'][zline,yline,100])
            str_data.n = f.variables['Density'][:,yline,xline]
            str_data.Vx = f.variables['Ux'][:,yline,xline]
            str_data.Vy = f.variables['Uy'][:,yline,xline]
            str_data.Vz = f.variables['Uz'][:,yline,xline] 
            str_data.T = f.variables['Temperature'][:,yline,xline]
            
    f.close()

def readNetcdf_grid(filepath,filename,str_data) :

    print(filepath+filename)
    f = nc.Dataset(filepath+filename,'r')
##    for v in f.variables:
##        print(v)


    str_data.V_sw = f.variables['vxs'][:]
    str_data.gstep = f.variables['gstep'][:]
    str_data.pos_planet = f.variables['s_centr'][:]
    str_data.c_omegapi = f.variables['phys_length'][:]
    str_data.s_min = f.variables['s_min'][:]
    str_data.s_max = f.variables['s_max'][:]
    
    str_data.x  = f.variables['X_axis'][:] 
    str_data.y  = f.variables['Y_axis'][:]
    str_data.z  = f.variables['Z_axis'][:]
    str_data.ref_dens = f.variables['phys_density'][:]

    f.close()

def read_file_particle(filepath,filename,str_data) :
  
    f = nc.Dataset(filepath+filename,'r')

##    for v in f.variables:
##        print(v)
    
    str_data.gstep = f.variables['gstep'][:]
    str_data.pos_planet = f.variables['s_centr'][:]
    str_data.c_omegapi = f.variables['phys_length'][:]
    str_data.s_min = f.variables['s_min'][:]
    
    str_data.x = -(f.variables['particule_x'][:]*str_data.gstep[0] - str_data.pos_planet[0]) 
    str_data.y = -(f.variables['particule_y'][:]*str_data.gstep[1] - str_data.pos_planet[1]) 
    str_data.z = f.variables['particule_z'][:]*str_data.gstep[2] - str_data.pos_planet[2] 

    str_data.Vx = -f.variables['particule_vx'][:] 
    str_data.Vy = -f.variables['particule_vy'][:] 
    str_data.Vz = f.variables['particule_vz'][:] 

    f.close()
    
