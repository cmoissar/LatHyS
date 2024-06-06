import pylab as pl
import scipy as scipy
import numpy as np
import matplotlib
import math
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.gridspec as gridspec

from class_data import *
from read_netcdf import *
from magfield_functions import *

from pylab import *

import pdb

matplotlib.rcParams['xtick.direction'] = 'out'
matplotlib.rcParams['ytick.direction'] = 'out'

matplotlib.rc('font',**{'family':'serif'})
plt.rc('text', usetex=True)
plt.rc('font')
matplotlib.rc('text', usetex=True)
matplotlib.rcParams['text.latex.preamble'] = [r'\boldmath']

fsize = 45

def xy(r,phi):
  return r*np.cos(phi), r*np.sin(phi)

def plot_color_maps(cl_data,data2plot,filepath_out,min_value=0.,max_value=0.,
                    sc_location = [-1,-1,-1],lang='En',flowline=0,xflow=0,
                    yflow=0,start_line=[0,0,0],end_line=[0,0,0]) :
    
    print(data2plot)
    
    pl.close()
    plt.cla()
        
    [XmeshXY,YmeshXY] = scipy.meshgrid(cl_data.x,cl_data.y)
    [YmeshYZ,ZmeshYZ] = scipy.meshgrid(cl_data.y,cl_data.z)    
    [XmeshXZ,ZmeshXZ] = scipy.meshgrid(cl_data.x,cl_data.z)  
    
    Xmin = int(XmeshXY.min())+2
    Xmax = int(XmeshXY.max())+2
    Ymin = int(YmeshXY.min())+2
    Ymax = int(YmeshXY.max())+2
    Zmin = int(ZmeshXZ.min())-1
    Zmax = int(ZmeshXZ.max())-1
    

    if data2plot == 'B':
        title0 = 'Magnetic field strength'
        data_name = 'Bfield'
        data_XY = cl_data.XY.B
        data_XZ = cl_data.XZ.B
        data_YZ = cl_data.YZ.B
        cb_title = 'B (nT)'
        cmap = plt.get_cmap('plasma')

    elif data2plot == 'Bx':
      title0 = 'Bx'
      data_name = 'Bx'
      data_XY = cl_data.XY.Bx
      data_XZ = cl_data.XZ.Bx
      data_YZ = cl_data.YZ.Bx

      cb_title = 'Bx (nT)'
      cmap = plt.get_cmap('seismic')
      
    elif data2plot == 'Bz':
      title0 = 'Bz'
      data_name = 'Bz'
      data_XY = cl_data.XY.Bz
      data_XZ = cl_data.XZ.Bz
      data_YZ = cl_data.YZ.Bz
      obstacleXY = np.where(np.sqrt(XmeshXY*XmeshXY + YmeshXY*YmeshXY) <= 20.)        
      data_XY[obstacleXY] = 0.

      obstacleXZ = np.where(np.sqrt(XmeshXZ*XmeshXZ + ZmeshXZ*ZmeshXZ) <= 20.)        
      data_XZ[obstacleXZ] = 0.

      obstacleYZ = np.where(np.sqrt(YmeshYZ*YmeshYZ + ZmeshYZ*ZmeshYZ) <= 20.)        
      data_YZ[obstacleYZ] = 0.
      cb_title = 'Bz (nT)'
      cmap = plt.get_cmap('seismic')
      
    elif data2plot == 'Pdyn':
        if(lang == 'En'):
          title0 = 'Dynamic pressure'
        elif(lang == 'Fr'):
          title0 = 'Pression dynamique'
        data_name = 'Pdyn'
        data_XY = cl_data.XY.Pdyn
        data_XZ = cl_data.XZ.Pdyn
        data_YZ = cl_data.YZ.Pdyn

        cb_title = '(nPa)'
        cmap = plt.get_cmap('plasma')
        
    elif data2plot == 'Pmag':
      if(lang == 'En'):
        title0 = 'Magnetic pressure'
      elif(lang == 'Fr'):
        title0 = r'Pression magn\'{e}tique'
      data_name = 'Pmag'
      data_XY = cl_data.XY.Pmag
      data_XZ = cl_data.XZ.Pmag
      data_YZ = cl_data.YZ.Pmag
      cb_title = '(nPa)'
      cmap = plt.get_cmap('plasma')
      
    elif data2plot == 'Pth':
      if(lang == 'En'):
        title0 = 'Thermal pressure'
      elif(lang == 'Fr'):
        title0 = 'Pression thermique'
      
      data_name = 'Pth'
      data_XY = cl_data.XY.Pth
      data_XZ = cl_data.XZ.Pth
      data_YZ = cl_data.YZ.Pth
      cb_title = '(nPa)'
      cmap = plt.get_cmap('plasma')
      
    elif data2plot == 'Pel':
        title0 = 'Pressure due to the electric field'
        data_name = 'Pel'
        data_XY = cl_data.XY.Pel
        data_XZ = cl_data.XZ.Pel
        data_YZ = cl_data.YZ.Pel
        cmap = plt.get_cmap('plasma')
        
    elif data2plot == 'Ptot':
        if(lang == 'En'):
          title0 = 'Total pressure'
        elif(lang == 'Fr'):
          title0 = 'Pression totale'
        data_name = 'Ptot'
        data_XY = cl_data.XY.Ptot
        data_XZ = cl_data.XZ.Ptot
        data_YZ = cl_data.YZ.Ptot
        cb_title = 'Ptot (nPa)'
        cmap = plt.get_cmap('plasma')
        
        
    elif data2plot == 'n':
      if lang == 'En':
        title0 = 'Ion density'
      elif lang == 'Fr':
        title0 = r"Densit\'{e} d'ions"
      
      data_name = 'dens'
      data_XY = cl_data.XY.n
      data_XZ = cl_data.XZ.n
      data_YZ = cl_data.YZ.n
      cb_title = r"\boldmath{$(\mathrm{cm}^{-3})$}"
      cmap = plt.get_cmap('plasma')
        
    elif data2plot == 'V':
      if lang == 'Fr':
        title0 = 'Vitesse des ions'
      else :
        title0 = 'Ion velocity'
        
      data_name = 'Vtot'
      data_XY = cl_data.XY.V
      data_XZ = cl_data.XZ.V
      data_YZ = cl_data.YZ.V
      cb_title = 'V (km/s)'
      cmap = plt.get_cmap('plasma')

    elif data2plot == 'Vx':
      if lang == 'En':
        title0 = 'Ion velocity along x'
      elif lang == 'Fr':
        title0 = 'Vitesse des ions (suivant x)'
              
      data_name = 'Vx'
      data_XY = cl_data.XY.Vx
      data_XZ = cl_data.XZ.Vx
      data_YZ = cl_data.YZ.Vx
      cb_title = 'Vx (km/s)'
      cmap = plt.get_cmap('seismic')

    elif data2plot == 'Tp':
      if lang == 'Fr':
        title0 = r'Temp\'{e}rature des ions'
      else:
        title0 = 'Ion temperature'
      data_name = 'Tp'
      data_XY = cl_data.XY.T*11605
      data_XZ = cl_data.XZ.T*11605
      data_YZ = cl_data.YZ.T*11605
      cb_title = '(K)'
      cmap = plt.get_cmap('plasma')
      
    elif data2plot == 'Tmag':
      if lang == 'En':
        title0 = 'Magnetic tension'
      elif lang == 'Fr':
        title0 = r'Tension magn\'{e}tique'
        
      data_name = 'Tmag'
      data_XY = cl_data.XY.Tmag
      data_XZ = cl_data.XZ.Tmag
      data_YZ = cl_data.YZ.Tmag
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('plasma')
      
    elif data2plot == 'Tmagx':
      if lang == 'En':
        title0 = 'Magnetic tension along x'
      elif lang == 'Fr':
        title0 = r'Tension magn\'{e}tique suivant x'
        
      data_name = 'Tmagx'
      data_XY = cl_data.XY.Tmagx
      data_XZ = cl_data.XZ.Tmagx
      data_YZ = cl_data.YZ.Tmagx
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')
        
    elif data2plot == 'Tmagy':
      if lang == 'En':
        title0 = 'Magnetic tension along y'
      elif lang == 'Fr':
        title0 = r'Tension magn\'{e}tique suivant y'
      data_name = 'Tmagy'
      data_XY = cl_data.XY.Tmagy
      data_XZ = cl_data.XZ.Tmagy
      data_YZ = cl_data.YZ.Tmagy
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')
        
    elif data2plot == 'Tmagz':
      if lang == 'En':
        title0 = 'Magnetic tension along z'
      elif lang == 'Fr':
        title0 = r'Tension magn\'{e}tique suivant z'
      data_name = 'Tmagz'
      data_XY = cl_data.XY.Tmagz
      data_XZ = cl_data.XZ.Tmagz
      data_YZ = cl_data.YZ.Tmagz
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')
        
    elif data2plot == 'dB2':
        title0 = 'Derivative of B2'
        data_name = 'dB2'
        data_XY = cl_data.XY.dB2
        data_XZ = cl_data.XZ.dB2
        data_YZ = cl_data.YZ.dB2
        cmap = plt.get_cmap('plasma')
        
    elif data2plot == 'Beta':
        title0 = 'Plasma beta'
        data_name = 'beta'
        data_XY = cl_data.XY.beta
        data_YZ = cl_data.YZ.beta
        data_XZ = cl_data.XZ.beta
        cmap = plt.get_cmap('plasma')
        cb_title = ''
        
    elif data2plot == 'Grad_Pmag':
      if lang == 'En':
        title0 = 'Magnetic pressure gradient'
      elif lang == 'Fr':
        title0 = r'Gradient de pression magn\'{e}tique'
        
      data_name = 'grad_Pmag'
      data_XY = cl_data.XY.grad_Pmag
      data_XZ = cl_data.XZ.grad_Pmag
      data_YZ = cl_data.YZ.grad_Pmag
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('plasma')
      
    elif data2plot == 'Grad_Pmagx':
      if lang == 'En':
        title0 = 'Magnetic pressure gradient along x'
      elif lang == 'Fr':
        title0 = r'\boldmath{$-\partial\left(B^2/2\mu_{0}\right)/\partial x $}'
  
      data_name = 'grad_Pmagx'
      data_XY = cl_data.XY.grad_Pmagx
      data_XZ = cl_data.XZ.grad_Pmagx
      data_YZ = cl_data.YZ.grad_Pmagx
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')
      
    elif data2plot == 'Grad_Pmagy':
      if lang == 'En':
        title0 = 'Magnetic pressure gradient along x'
      elif lang == 'Fr':
        title0 = r'Gradient de pression magn\'{e}tique (suivant -y)'
        
      data_name = 'grad_Pmagy'
      data_XY = cl_data.XY.grad_Pmagy
      data_XZ = cl_data.XZ.grad_Pmagy
      data_YZ = cl_data.YZ.grad_Pmagy
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')
      
    elif data2plot == 'Grad_Pmagz':
      if lang == 'En':
        title0 = 'Magnetic pressure gradient along x'
      elif lang == 'Fr':
        title0 = r'Gradient de pression magn\'{e}tique (suivant -z)'
        
      data_name = 'grad_Pmagz'
      data_XY = cl_data.XY.grad_Pmagz
      data_XZ = cl_data.XZ.grad_Pmagz
      data_YZ = cl_data.YZ.grad_Pmagz
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')

    elif data2plot == 'Grad_Pth':
      if lang == 'En':
        title0 = 'Thermal pressure gradient'
      elif lang == 'Fr':
        title0 = 'Gradient de pression thermique'
        
      data_name = 'grad_Pth'
      data_XY = cl_data.XY.grad_Pth
      data_XZ = cl_data.XZ.grad_Pth
      data_YZ = cl_data.YZ.grad_Pth
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('plasma')

    elif data2plot == 'Grad_Pthx':
      if lang == 'En':
        title0 = 'Gradient along x of the thermal pressure'
      elif lang == 'Fr':
        title0 = r'\boldmath{$- \partial p/\partial x$}'
        
      data_name = 'grad_Pthx'
      data_XY = cl_data.XY.grad_Pthx
      data_XZ = cl_data.XZ.grad_Pthx
      data_YZ = cl_data.YZ.grad_Pthx
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')

    elif data2plot == 'Grad_Pthy':
      if lang == 'En':
        title0 = 'Gradient along y of the thermal pressure'
      elif lang == 'Fr':
        title0 = 'Gradient de pression thermique (suivant y)'
        
      data_name = 'grad_Pthy'
      data_XY = cl_data.XY.grad_Pthy
      data_XZ = cl_data.XZ.grad_Pthy
      data_YZ = cl_data.YZ.grad_Pthy
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')

    elif data2plot == 'Grad_Pthz':
      if lang == 'En':
        title0 = 'Gradient along z of the thermal pressure'
      elif lang == 'Fr':
        title0 = 'Gradient de pression thermique (suivant z)'
        
      data_name = 'grad_Pthz'
      data_XY = cl_data.XY.grad_Pthz
      data_XZ = cl_data.XZ.grad_Pthz
      data_YZ = cl_data.YZ.grad_Pthz
      cb_title = r"\boldmath{$(\mathrm{nPa}.\mathrm{km}^{-1})$}"
      cmap = plt.get_cmap('seismic')
      
    elif data2plot == 'psi':
        if(lang == 'En'):
          title0 = r"\boldmath{$\psi$}"
        elif(lang == 'Fr'):
          title0 = r"\boldmath{$\psi$}"
        data_name = 'psi'
        data_XY = cl_data.XY.psi
        data_XZ = cl_data.XZ.psi
        data_YZ = cl_data.YZ.psi
        cb_title = r"\boldmath{$(^{\circ})}$}"
        cmap = plt.get_cmap('plasma')
        
    elif data2plot == 'Etot':
        title0 = 'Electric field'
        data_name = 'Etot'
        data_XY = cl_data.XY.E
        data_XZ = cl_data.XZ.E
        data_YZ = cl_data.YZ.E
        cmap = plt.get_cmap('plasma')

    elif data2plot == 'Ma':
      if lang == 'Fr':
        title0 = 'Nombre de Mach d''Alfven'
      else :
        title0 = 'Alfven Mach number'
        
      data_name = 'Ma'
      data_XY = cl_data.XY.Ma
      data_XZ = cl_data.XZ.Ma
      data_YZ = cl_data.YZ.Ma
      cb_title = ''
      cmap = plt.get_cmap('plasma')

    elif data2plot == 'Ms':
      if lang == 'Fr':
        title0 = 'Nombre de Mach sonique'
      else :
        title0 = 'Sonic Mach number'
        
      data_name = 'Ms'
      data_XY = cl_data.XY.Ms
      data_XZ = cl_data.XZ.Ms
      data_YZ = cl_data.YZ.Ms
      cb_title = ''
      cmap = plt.get_cmap('plasma')

    elif data2plot == 'Mms':
      if lang == 'Fr':
        title0 = r'Nombre de Mach magn\'{e}tosonore'
      else :
        title0 = 'Magnetosonic Mach number'
        
      data_name = 'Mms'
      data_XY = cl_data.XY.Mms
      data_XZ = cl_data.XZ.Mms
      data_YZ = cl_data.YZ.Mms
      cb_title = ''
      cmap = plt.get_cmap('plasma')
      
    if min_value == 0 and max_value == 0:
        datamin = np.array([np.nanmin(data_XY),np.nanmin(data_XZ),np.nanmin(data_YZ)])
        datamax = np.array([np.nanmax(data_XY),np.nanmax(data_XZ),np.nanmax(data_YZ)])    
        min_value = int(datamin.min())
        max_value = int(datamax.max())

        #min_value = int(np.nanmean(data_XY)+np.nanmean(data_XZ)+np.nanmean(data_YZ)) / 6.0
        max_value = int(np.nanmean(data_XY)+np.nanmean(data_XZ)+np.nanmean(data_YZ))
        # n = max_value
        # max_value = int(n / 10 ** (int(math.log10(n)))) * 10 ** (int(math.log10(n)))

    # Number of color levels
    levels = MaxNLocator(nbins=255).tick_values(min_value, max_value)

    # Chosen colormap
    #cmap = plt.get_cmap('jet') #'jet'

    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

#    gs = gridspec.GridSpec(1, 3, width_ratios=[1,
#        float(cl_data.y.shape[0])/cl_data.x.shape[0],1],left=0.05,right=0.9)#,wspace=0.6) 
    
    gs = gridspec.GridSpec(1,3)

    fig = pl.figure(figsize=[25,13.2])

    plt.suptitle(title0 +' (t = '+cl_data.time+
              r"\boldmath{$\, \Omega_{\mathrm{ci}}^{-1}$}"+')',
                 fontsize=fsize*1.3,y=0.95,weight="bold")

    sub1 = plt.subplot(gs[0])
    fig = pl.pcolor(XmeshXY,YmeshXY,data_XY[:,:],cmap=cmap,norm=norm)

    plt.subplot(sub1).set_aspect('equal')
    plt.axis([cl_data.x[0],cl_data.x[-1],
              cl_data.y[-1],cl_data.y[0]])

    plt.xlabel('x',fontsize=24)
    plt.ylabel('y',fontsize=24)
    #plot only the magnetosheath
    plt.xlim([-2*cl_data.x[-1],cl_data.x[-1]])
    plt.ylim([cl_data.y[-1],cl_data.y[0]])

    plt.setp(sub1.get_xticklabels(),fontsize=18)
    plt.setp(sub1.get_yticklabels(),fontsize=18)

    tick_locs = MaxNLocator(nbins=5).tick_values(cl_data.x[-1], -2*cl_data.x[-1])
    plt.xticks(tick_locs, [r"$\mathbf{%s}$" % int(x) for x in tick_locs])
    tick_locs = MaxNLocator(nbins=5).tick_values(min(cl_data.y), max(cl_data.y))
    plt.yticks(tick_locs, [r"$\mathbf{%s}$" % int(x) for x in tick_locs])

    if flowline == 1:
      sub1.plot(xflow,yflow,c='w',ls='-' ,linewidth=2.5)

    if not np.all([start_line==end_line]):
      sub1.plot([start_line[0],end_line[0]],[start_line[1],end_line[1]],
                c='w',ls='-' ,linewidth=2.5)

    sub2 = plt.subplot(gs[1])
    fig = pl.pcolor(YmeshYZ,ZmeshYZ,data_YZ[:,:],cmap=cmap,norm=norm)

    plt.subplot(sub2).set_aspect('equal')    
    plt.axis([cl_data.y[0],cl_data.y[cl_data.y.shape[0]-1],
        cl_data.z[0],cl_data.z[cl_data.z.shape[0]-1]])

    plt.xlabel('y',fontsize=24)
    plt.ylabel('z',fontsize=24)
    plt.xlim([cl_data.y[cl_data.y.shape[0]-1],cl_data.y[0]])
    plt.ylim([cl_data.z[0],cl_data.z[cl_data.z.shape[0]-1]])

    plt.setp(sub2.get_xticklabels(),fontsize=18)
    plt.setp(sub2.get_yticklabels(),fontsize=18)


    tick_locs = MaxNLocator(nbins=5).tick_values(min(cl_data.y), max(cl_data.y))
    plt.yticks(tick_locs, [r"$\mathbf{%s}$" % int(x) for x in tick_locs])
    tick_locs = MaxNLocator(nbins=5).tick_values(min(cl_data.z), max(cl_data.z))
    plt.yticks(tick_locs, [r"$\mathbf{%s}$" % int(x) for x in tick_locs])

    if not np.all([start_line==end_line]):
      sub2.plot([start_line[1],end_line[1]],[start_line[2],end_line[2]],
                c='w',ls='-' ,linewidth=2.5)
    
    sub3 = plt.subplot(gs[2])
    fig = pl.pcolor(XmeshXZ,ZmeshXZ,data_XZ[:,:],cmap=cmap,norm=norm)

    #plt.subplot(sub1).set_aspect(aspect=4) #Used for long box
    plt.subplot(sub3).set_aspect('equal')
    plt.axis([cl_data.x[0],cl_data.x[-1],
        cl_data.z[-1],cl_data.z[0]])

    plt.xlabel('x',fontsize=24)
    plt.ylabel('z',fontsize=24)
    plt.xlim([-2*cl_data.x[-1],cl_data.x[-1]])
    plt.ylim([cl_data.y[-1],cl_data.y[0]])

    plt.setp(sub3.get_xticklabels(),fontsize=18)
    plt.setp(sub3.get_yticklabels(),fontsize=18)

    tick_locs = MaxNLocator(nbins=5).tick_values(cl_data.x[-1], -2*cl_data.x[-1])
    plt.xticks(tick_locs, [r"$\mathbf{%s}$" % int(x) for x in tick_locs])
    tick_locs = MaxNLocator(nbins=5).tick_values(min(cl_data.z), max(cl_data.z))
    plt.yticks(tick_locs, [r"$\mathbf{%s}$" % int(x) for x in tick_locs])

    if not np.all([start_line==end_line]):
      sub3.plot([start_line[0],end_line[0]],[start_line[2],end_line[2]],
                c='w',ls='-' ,linewidth=2.5)
      

#
#    plt.subplot(gs[2:])    
#    #moyenne arithmétique un peu bête. Ne prend pas en compte la taille en Y ou Z
#    signal = (np.average(data_XY[:-1,:],0)+np.average(data_XZ[:-1,:],0)) / 2
#    plt.plot(cl_data.x, signal )
#    plt.gca().invert_xaxis()
#    plt.xlabel('x',fontsize=fsize)
#    tick_locs = range(Xmin, Xmax, int((Xmax-Xmin)/5))
#    plt.xticks(tick_locs, [r"$\mathbf{%s}$" % x for x in tick_locs])

    cax = pl.axes([0.92, 0.138, 0.025, 0.725])

    if data2plot == 'B' or data2plot == 'Bx' :
      plt.text(-0.15, 1.02, cb_title, fontsize=fsize/1.5)
      tick_locs = np.linspace(min_value,max_value,5)
      cb = plt.colorbar(cax=cax,ticks=tick_locs)
    elif data2plot == 'V':
      tick_locs = np.linspace(min_value,max_value,5)
      cb = plt.colorbar(cax=cax,ticks=tick_locs)
    elif data2plot == 'n':
      tick_locs = np.linspace(min_value,max_value,5)
      cb = plt.colorbar(cax=cax,ticks=tick_locs)
    elif data2plot == 'psi':
      plt.text(0.18, 1.02, cb_title, fontsize=fsize,weight='bold')
      if max_value == 60:
        tick_locs = [i for i in range(min_value,max_value+10,10)] 
        cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                       num=7),format='%3i')
      elif max_value == 90:
        tick_locs = [i for i in range(min_value,max_value+10,10)] 
        cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                     num=10),format='%3i')
      elif max_value == 100:
        tick_locs = [i for i in range(min_value,max_value+10,20)] 
        cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                     num=6),format='%3i')
      else:
        tick_locs = [i for i in range(min_value,max_value+10,10)] 
        cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                     num=int(max_value/10)+1))
   

    elif data2plot == 'Tmagx' or data2plot == 'Tmagy' or data2plot == 'Tmagz' or\
     data2plot == 'Grad_Pmagx' or data2plot == 'Grad_Pmagy' or data2plot == 'Grad_Pmagz' or\
     data2plot == 'Grad_Pthx' or data2plot == 'Grad_Pthy' or data2plot == 'Grad_Pthz':
      tick_locs = [i for i in range(int(min_value*1.e+5),int(max_value*1.e+5)+8,8)]
      
      cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                   num=11),format='%8.1e')
      plt.text(-0.5, 1.02, cb_title, fontsize=fsize,weight='bold')
      plt.text(1.5, 0.96, r'\boldmath{$\times 10^{-5}$}', fontsize=fsize,weight='bold')

    elif data2plot == 'Tmag' or data2plot == 'Grad_Pmag' or data2plot == 'Grad_Pth':
      tick_locs = np.linspace(min_value*1.e+5,max_value*1.e+5,num=9)
      
      cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                   num=9),format='%8.1e')
      plt.text(-0.5, 1.02, cb_title, fontsize=fsize,weight='bold')
      plt.text(1.5, 0.96, r'\boldmath{$\times 10^{-5}$}', fontsize=fsize,weight='bold')
      

    elif data2plot == 'Pmag' or data2plot == 'Pth' or data2plot == 'Pdyn' or data2plot == 'Ptot':
      plt.text(-0.25, 1.02, cb_title, fontsize=fsize,weight='bold')

      tick_locs = np.linspace(min_value,max_value,num=6,endpoint=True)
      cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                   num=6,endpoint=True))
      
    elif data2plot == 'Tp':
      tick_locs = np.linspace(min_value,max_value,num=11)
      cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                   num=11))
      plt.text(0, 1.02, cb_title, fontsize=fsize,weight='bold')     

    else:
      cb = plt.colorbar(cax=cax,ticks=np.linspace(min_value,max_value,
                                                  num=11))
      plt.text(-0.25, 1.02, cb_title, fontsize=fsize,weight='bold')
      tick_locs = np.linspace(min_value,max_value,num=11)

#    cb.ax.set_yticklabels(tick_locs, [r"$\mathbf{%s}$" % x for x in tick_locs])

    for t in cb.ax.get_yticklabels():
       t.set_fontsize(fsize)

    if not np.all([start_line==end_line]):

      fileout = data_name+'_trajectory('+str(int(start_line[0])) + \
              '_'+str(int(start_line[1]))+'_'+str(int(start_line[2])) \
              +')_('+str(int(end_line[0])) + '_'+str(int(end_line[1]))+ \
              '_'+str(int(end_line[2])) + ')_t='+cl_data.time+'.png'
      
      print(filepath_out+fileout)
      plt.savefig(filepath_out+fileout,dpi=300)
      print('Figure saved as ',filepath_out+fileout)
    else :
      print(filepath_out+data_name+'_t='+cl_data.time+'.png')
      plt.savefig(filepath_out+data_name+'_t='+cl_data.time+'.png',dpi=300)
      print('Figure saved as ',filepath_out+data_name+'_t='+cl_data.time+'.png')

    return
