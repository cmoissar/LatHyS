import numpy as np
import math as m

from class_data import *
from energy_density import *

def calculate_norm(magfield):

    ## Calculates the magnitude of the magnetic field
    magfield.B = np.sqrt(np.add(np.add(np.multiply(magfield.Bx,magfield.Bx),
        np.multiply(magfield.By,magfield.By)),np.multiply(magfield.Bz,magfield.Bz)))
        
def calculate_theta(magfield):
    ## Calculates the angle between the magnetic field vector and the z direction
    magfield.theta = np.arccos(np.divide(magfield.Bz,magfield.B))
        
def calculate_phi(magfield):
   ## Calculates the angle between the projection of the magnetic field vector
   ## in the XY plane and the x axis.
   sign_By = np.ones(magfield.Bx.shape)
   sign_By = np.copysign(sign_By,magfield.By)

   magfield.phi = np.multiply(np.arccos(np.divide(magfield.Bx,np.sqrt(np.add(np.multiply(magfield.Bx,magfield.Bx),
    np.multiply(magfield.By,magfield.By))))),sign_By) 

def calculate_clockangle(magfield):
   ## Calculates the angle between the projection of the magnetic field vector
   ## in the YZ plane and the z axis.

   sign_By = np.ones(magfield.Bx.shape)
   sign_By = np.copysign(sign_By,magfield.By)
   
   magfield.clockangle = np.multiply(np.arccos(np.divide(magfield.Bz,np.sqrt(np.add(np.multiply(magfield.By,magfield.By),
   np.multiply(magfield.Bz,magfield.Bz))))),sign_By) 

def calculate_coneangle(magfield):
    ## Calculates the angle between the magnetic field vector and
    ## the x axis.

    magfield.coneangle = np.arccos(np.divide(magfield.Bx,magfield.B))
    

def calculate_mag_tension(magfield,magfieldsup,magfieldinf):

    ## Calculates the magnetic tension in the three considered planes of the simulation
    ## Values of the magnetic tension given in ??
    
    print('Inside the calc_mag_tension function')
    magfield.XY.Tmag  = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.Tmagx = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.Tmagy = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.Tmagz = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))

    magfield.XZ.Tmag  = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.Tmagx = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.Tmagy = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.Tmagz = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    
    magfield.YZ.Tmag  = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.Tmagx = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.Tmagy = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.Tmagz = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))

    mu_0 = np.pi*4.e-7
    
    print('First loop')
    for i in range(1,magfield.x.shape[0]-1):
        for j in range(1,magfield.y.shape[0]-1):
            magfield.XY.Tmagx[j,i] = \
                 (magfield.XY.Bx[j,i]*(magfield.XY.Bx[j,i-1]-magfield.XY.Bx[j,i+1])/ \
                 (magfield.gstep[2]*magfield.c_omegapi)) + \
                 (magfield.XY.By[j,i]*(magfield.XY.Bx[j-1,i]-magfield.XY.Bx[j+1,i])/ \
                 (magfield.gstep[1]*magfield.c_omegapi)) + \
                 (magfield.XY.Bz[j,i]*(magfieldsup.XY.Bx[j,i]-magfieldinf.XY.Bx[j,i])/ \
                 (magfield.gstep[0]*magfield.c_omegapi))

            magfield.XY.Tmagy[j,i] = \
                (magfield.XY.Bx[j,i]*(magfield.XY.By[j,i-1]-magfield.XY.By[j,i+1])/ \
                (magfield.gstep[2]*magfield.c_omegapi)) + \
                (magfield.XY.By[j,i]*(magfield.XY.By[j-1,i]-magfield.XY.By[j+1,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi)) + \
                (magfield.XY.Bz[j,i]*(magfieldsup.XY.By[j,i]-magfieldinf.XY.By[j,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))

            magfield.XY.Tmagz[j,i] = \
                 (magfield.XY.Bx[j,i]*(magfield.XY.Bz[j,i-1]-magfield.XY.Bz[j,i+1])/ \
                  (magfield.gstep[2]*magfield.c_omegapi)) + \
                 (magfield.XY.By[j,i]*(magfield.XY.Bz[j-1,i]-magfield.XY.Bz[j+1,i])/ \
                  (magfield.gstep[1]*magfield.c_omegapi)) + \
                  (magfield.XY.Bz[j,i]*(magfieldsup.XY.Bz[j,i]-magfieldinf.XY.Bz[j,i])/ \
                  (magfield.gstep[0]*magfield.c_omegapi))       

            magfield.XY.Tmagx[j,i] = magfield.XY.Tmagx[j,i]*1.e-9/(2.*mu_0)
            magfield.XY.Tmagy[j,i] = magfield.XY.Tmagy[j,i]*1.e-9/(2.*mu_0)
            magfield.XY.Tmagz[j,i] = magfield.XY.Tmagz[j,i]*1.e-9/(2.*mu_0)
                                                             
            magfield.XY.Tmag[j,i] = np.sqrt(magfield.XY.Tmagx[j,i]*magfield.XY.Tmagx[j,i] \
                                           + magfield.XY.Tmagy[j,i]*magfield.XY.Tmagy[j,i] \
                                           + magfield.XY.Tmagz[j,i]*magfield.XY.Tmagz[j,i])  
    print('First loop end - starting second loop')
    #print(magfield.XY.Tmag[0:20,50])
    #print(magfield.c_omegapi)
    #print(magfield.gstep)
    
    for i in range(1,magfield.x.shape[0]-1):
        for j in range(1,magfield.z.shape[0]-1):
            magfield.XZ.Tmagx[j,i] = \
                (magfield.XZ.Bx[j,i]*(magfield.XZ.Bx[j,i-1]-magfield.XZ.Bx[j,i+1])/ \
                (magfield.gstep[2]*magfield.c_omegapi)) + \
                (magfield.XZ.By[j,i]*(magfieldsup.XZ.Bx[j,i]-magfieldinf.XZ.Bx[j,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi)) + \
                (magfield.XZ.Bz[j,i]*(magfield.XZ.Bx[j+1,i]-magfield.XZ.Bx[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))

            magfield.XZ.Tmagy[j,i] = \
                (magfield.XZ.Bx[j,i]*(magfield.XZ.By[j,i-1]-magfield.XZ.By[j,i+1])/ \
                (magfield.gstep[2]*magfield.c_omegapi)) + \
                (magfield.XZ.By[j,i]*(magfieldsup.XZ.By[j,i]-magfieldinf.XZ.By[j,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi)) + \
                (magfield.XZ.Bz[j,i]*(magfield.XZ.By[j+1,i]-magfield.XZ.By[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))

            magfield.XZ.Tmagz[j,i] = \
                (magfield.XZ.Bx[j,i]*(magfield.XZ.Bz[j,i-1]-magfield.XZ.Bz[j,i+1])/ \
                (magfield.gstep[2]*magfield.c_omegapi)) + \
                (magfield.XZ.By[j,i]*(magfieldsup.XZ.Bz[j,i]-magfieldinf.XZ.Bz[j,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi)) + \
                (magfield.XZ.Bz[j,i]*(magfield.XZ.Bz[j+1,i]-magfield.XZ.Bz[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))

            magfield.XZ.Tmagx[j,i] = magfield.XZ.Tmagx[j,i]*1.e-9/(2.*mu_0)
            magfield.XZ.Tmagy[j,i] = magfield.XZ.Tmagy[j,i]*1.e-9/(2.*mu_0)
            magfield.XZ.Tmagz[j,i] = magfield.XZ.Tmagz[j,i]*1.e-9/(2.*mu_0)
                                                             
            magfield.XZ.Tmag[j,i] = np.sqrt(magfield.XZ.Tmagx[j,i]*magfield.XZ.Tmagx[j,i] \
                                           + magfield.XZ.Tmagy[j,i]*magfield.XZ.Tmagy[j,i] \
                                           + magfield.XZ.Tmagz[j,i]*magfield.XZ.Tmagz[j,i])  
    print('Second loop end - starting third loop')        
    #print(magfield.XZ.Tmag[0:20,50])

    for i in range(1,magfield.y.shape[0]-1):
        for j in range(1,magfield.z.shape[0]-1):
            magfield.YZ.Tmagx[j,i] = \
                (magfield.YZ.Bx[j,i]*(magfieldsup.YZ.Bx[j,i]-magfield.YZ.Bx[j,i])/ \
                (magfield.gstep[2]*magfield.c_omegapi)) + \
                (magfield.YZ.By[j,i]*(magfield.YZ.Bx[j,i-1]-magfield.YZ.Bx[j,i+1])/ \
                (magfield.gstep[1]*magfield.c_omegapi)) + \
                (magfield.YZ.Bz[j,i]*(magfield.YZ.Bx[j+1,i]-magfield.YZ.Bx[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))

            magfield.YZ.Tmagy[j,i] = \
                (magfield.YZ.Bx[j,i]*(magfieldsup.YZ.By[j,i]-magfieldinf.YZ.By[j,i])/ \
                (magfield.gstep[2]*magfield.c_omegapi)) + \
                (magfield.YZ.By[j,i]*(magfield.YZ.By[j,i-1]-magfield.YZ.By[j,i+1])/ \
                (magfield.gstep[1]*magfield.c_omegapi)) + \
                (magfield.YZ.Bz[j,i]*(magfield.YZ.By[j+1,i]-magfield.YZ.By[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))

            magfield.YZ.Tmagz[j,i] = \
                (magfield.YZ.Bx[j,i]*(magfieldsup.YZ.Bz[j,i]-magfieldinf.YZ.Bz[j,i])/ \
                (magfield.gstep[2]*magfield.c_omegapi)) + \
                (magfield.YZ.By[j,i]*(magfield.YZ.Bz[j,i-1]-magfield.YZ.Bz[j,i+1])/ \
                (magfield.gstep[1]*magfield.c_omegapi)) + \
                (magfield.YZ.Bz[j,i]*(magfield.YZ.Bz[j+1,i]-magfield.YZ.Bz[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))      

            magfield.YZ.Tmagx[j,i] = magfield.YZ.Tmagx[j,i]*1.e-9/(2.*mu_0)
            magfield.YZ.Tmagy[j,i] = magfield.YZ.Tmagy[j,i]*1.e-9/(2.*mu_0)
            magfield.YZ.Tmagz[j,i] = magfield.YZ.Tmagz[j,i]*1.e-9/(2.*mu_0)
            
            magfield.YZ.Tmag[j,i] = np.sqrt(magfield.YZ.Tmagx[j,i]*magfield.YZ.Tmagx[j,i] \
                                           + magfield.YZ.Tmagy[j,i]*magfield.YZ.Tmagy[j,i] \
                                           + magfield.YZ.Tmagz[j,i]*magfield.YZ.Tmagz[j,i])  
         
    #print(magfield.YZ.Tmag[0:20,50])


    
def calculate_psi(magfield_SW,magfield_Msheath):

# Calculates the angle between the magnetic field direction in the solar wind
# and in the magnetosheath

    magfield_Msheath.psi = np.arccos(np.divide(np.add(np.add(
        np.multiply(magfield_SW.Bx,magfield_Msheath.Bx),
        np.multiply(magfield_SW.By,magfield_Msheath.By)),
        np.multiply(magfield_SW.Bz,magfield_Msheath.Bz)),
        np.multiply(magfield_SW.B,magfield_Msheath.B)))

def calculate_psi_v2(magfield_SW,magfield_Msheath,plane='XY'):


##    zero = np.where(magfield_Msheath.B < 1e-3)
##    magfield_Msheath.B[zero] = magfield_SW.B

##    whereNaN = np.isnan(magfield_Msheath.B)
##    magfield_Msheath.B[whereNaN] = magfield_SW.B
##
##    where_zero = np.where(np.absolute(magfield_Msheath.Bx) < 1.e-5 and
##                          np.absolute(magfield_Msheath.By) < 1.e-5 and
##                          np.absolute(magfield_Msheath.Bz) < 1.e-5)
##    
##    print(magfield_Msheath.Bx,magfield_Msheath.By,magfield_Msheath.Bz)
##    print(zero)
##    print(magfield_SW.B)
##    print(magfield_Msheath.B[220])

   # magfield_Msheath.psi[where_zero] = 0.
##    magfield_Msheath.psi = np.arccos(np.divide(np.add(np.add(
##        magfield_SW.Bx*magfield_Msheath.Bx,
##        magfield_SW.By*magfield_Msheath.By),
##        magfield_SW.Bz*magfield_Msheath.Bz),
##        magfield_SW.B*magfield_Msheath.B))                          
##    magfield_Msheath.psi = np.arccos(np.divide(np.add(np.add(
##        magfield_SW.Bx*magfield_Msheath.Bx,
##        magfield_SW.By*magfield_Msheath.By),
##        magfield_SW.Bz,magfield_Msheath.Bz),
##        magfield_SW.B*magfield_Msheath.B))

    if plane == 'XY':
        nx = magfield_Msheath.x.shape[0]
        ny = magfield_Msheath.y.shape[0]
    elif plane == 'XZ':
        nx = magfield_Msheath.x.shape[0]
        ny = magfield_Msheath.z.shape[0]
    elif plane == 'YZ':
        nx = magfield_Msheath.y.shape[0]
        ny = magfield_Msheath.z.shape[0]
        
    magfield_Msheath.XY.psi = np.zeros((magfield_Msheath.y.shape[0],
                                        magfield_Msheath.x.shape[0]))
    magfield_Msheath.YZ.psi = np.zeros((magfield_Msheath.z.shape[0],
                                        magfield_Msheath.y.shape[0]))
    magfield_Msheath.XZ.psi = np.zeros((magfield_Msheath.z.shape[0],
                                        magfield_Msheath.x.shape[0]))    

    for i in range(1,magfield_Msheath.x.shape[0]-1):
        for j in range(1,magfield_Msheath.y.shape[0]-1):
            magfield_Msheath.XY.psi[j,i] = np.arccos((
                magfield_SW.Bx*magfield_Msheath.XY.Bx[j,i] +
                magfield_SW.By*magfield_Msheath.XY.By[j,i] +
                magfield_SW.Bz*magfield_Msheath.XY.Bz[j,i])/(
                    magfield_SW.B*magfield_Msheath.XY.B[j,i]))

            #print(magfield_Msheath.XY.psi[j,i])

    for i in range(1,magfield_Msheath.y.shape[0]-1):
        for j in range(1,magfield_Msheath.z.shape[0]-1):
            magfield_Msheath.YZ.psi[j,i] = np.arccos((
                magfield_SW.Bx*magfield_Msheath.YZ.Bx[j,i] +
                magfield_SW.By*magfield_Msheath.YZ.By[j,i] +
                magfield_SW.Bz*magfield_Msheath.YZ.Bz[j,i])/(
                    magfield_SW.B*magfield_Msheath.YZ.B[j,i]))

            #(magfield_Msheath.YZ.psi[j,i])

    for i in range(1,magfield_Msheath.x.shape[0]-1):
        for j in range(1,magfield_Msheath.z.shape[0]-1):
            magfield_Msheath.XZ.psi[j,i] = np.arccos((
                magfield_SW.Bx*magfield_Msheath.XZ.Bx[j,i] +
                magfield_SW.By*magfield_Msheath.XZ.By[j,i] +
                magfield_SW.Bz*magfield_Msheath.XZ.Bz[j,i])/(
                    magfield_SW.B*magfield_Msheath.XZ.B[j,i]))

            #print(magfield_Msheath.XZ.psi[j,i])            
    

def calculate_mag_pressure_gradient(magfield,magfieldsup,magfieldinf):

    magfield.XY.grad_Pmag   = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.grad_Pmagx  = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.grad_Pmagy  = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.grad_Pmagz  = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    
    magfield.XZ.grad_Pmag   = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.grad_Pmagx  = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.grad_Pmagy  = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.grad_Pmagz  = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
        
    magfield.YZ.grad_Pmag   = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.grad_Pmagx  = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.grad_Pmagy  = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.grad_Pmagz  = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))    

    calc_mag_pressure(magfield.XY,calc_grad=1)
    calc_mag_pressure(magfield.XZ,calc_grad=1)
    calc_mag_pressure(magfield.YZ,calc_grad=1)
    
    calc_mag_pressure(magfieldsup.XY,calc_grad=1)
    calc_mag_pressure(magfieldsup.XZ,calc_grad=1)
    calc_mag_pressure(magfieldsup.YZ,calc_grad=1)    

    calc_mag_pressure(magfieldinf.XY,calc_grad=1)
    calc_mag_pressure(magfieldinf.XZ,calc_grad=1)
    calc_mag_pressure(magfieldinf.YZ,calc_grad=1)

##    magfield.XY.Pmag = magfield.XY.Pmag*1.e9
##    magfield.YZ.Pmag = magfield.YZ.Pmag*1.e9
##    magfield.XZ.Pmag = magfield.XZ.Pmag*1.e9
##
##    magfieldsup.XY.Pmag = magfieldsup.XY.Pmag*1.e9
##    magfieldsup.YZ.Pmag = magfieldsup.YZ.Pmag*1.e9
##    magfieldsup.XZ.Pmag = magfieldsup.XZ.Pmag*1.e9
##
##    magfieldinf.XY.Pmag = magfieldinf.XY.Pmag*1.e9
##    magfieldinf.YZ.Pmag = magfieldinf.YZ.Pmag*1.e9
##    magfieldinf.XZ.Pmag = magfieldinf.XZ.Pmag*1.e9
    
    for i in range(1,magfield.x.shape[0]-1):
        for j in range(1,magfield.y.shape[0]-1):

            magfield.XY.grad_Pmagx[j,i] = 1.e-9*(magfield.XY.Pmag[j,i-1] - magfield.XY.Pmag[j,i+1])/ \
                (magfield.gstep[2]*magfield.c_omegapi*2.)
            magfield.XY.grad_Pmagy[j,i] = 1.e-9*(magfield.XY.Pmag[j-1,i] - magfield.XY.Pmag[j+1,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi*2.)
            magfield.XY.grad_Pmagz[j,i] = 1.e-9*(magfieldsup.XY.Pmag[j,i] - magfieldinf.XY.Pmag[j,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi*2.)

            magfield.XY.grad_Pmag[j,i] = np.sqrt(magfield.XY.grad_Pmagx[j,i]*magfield.XY.grad_Pmagx[j,i] +
                                                 magfield.XY.grad_Pmagy[j,i]*magfield.XY.grad_Pmagy[j,i] +
                                                 magfield.XY.grad_Pmagz[j,i]*magfield.XY.grad_Pmagz[j,i])

    for i in range(1,magfield.x.shape[0]-1):
        for j in range(1,magfield.z.shape[0]-1):

            magfield.XZ.grad_Pmagx[j,i] = 1.e-9*(magfield.XZ.Pmag[j,i-1] - magfield.XZ.Pmag[j,i+1])/ \
                (magfield.gstep[2]*magfield.c_omegapi*2.)
            magfield.XZ.grad_Pmagy[j,i] = 1.e-9*(magfieldsup.XZ.Pmag[j,i] - magfieldinf.XZ.Pmag[j,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi*2.)
            magfield.XZ.grad_Pmagz[j,i] = 1.e-9*(magfield.XZ.Pmag[j+1,i] - magfield.XZ.Pmag[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi*2.)

            magfield.XZ.grad_Pmag[j,i] = np.sqrt(magfield.XZ.grad_Pmagx[j,i]*magfield.XZ.grad_Pmagx[j,i] +
                                                 magfield.XZ.grad_Pmagy[j,i]*magfield.XZ.grad_Pmagy[j,i] +
                                                 magfield.XZ.grad_Pmagz[j,i]*magfield.XZ.grad_Pmagz[j,i])

    for i in range(1,magfield.y.shape[0]-1):
        for j in range(1,magfield.z.shape[0]-1):

            magfield.YZ.grad_Pmagx[j,i] = 1.e-9*(magfieldsup.YZ.Pmag[j,i] - magfieldinf.YZ.Pmag[j,i])/ \
                (magfield.gstep[2]*magfield.c_omegapi*2.)
            magfield.YZ.grad_Pmagy[j,i] = 1.e-9*(magfield.YZ.Pmag[j,i-1] - magfield.YZ.Pmag[j,i+1])/ \
                (magfield.gstep[1]*magfield.c_omegapi*2.)
            magfield.YZ.grad_Pmagz[j,i] = 1.e-9*(magfield.YZ.Pmag[j+1,i] - magfield.YZ.Pmag[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi*2.)

            magfield.YZ.grad_Pmag[j,i] = np.sqrt(magfield.YZ.grad_Pmagx[j,i]*magfield.YZ.grad_Pmagx[j,i] +
                                                 magfield.YZ.grad_Pmagy[j,i]*magfield.YZ.grad_Pmagy[j,i] +
                                                 magfield.YZ.grad_Pmagz[j,i]*magfield.YZ.grad_Pmagz[j,i])

def calculate_dB2(magfield,magfieldsup,magfieldinf):

    ## Calculates the magnetic tension in the three considered planes of the simulation
    ## Values of the magnetic tension given in ????
    
    print('Inside the calc_mag_tension function')
    magfield.XY.dB2  = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.dB2x = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.dB2y = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))
    magfield.XY.dB2z = np.zeros((magfield.y.shape[0],magfield.x.shape[0]))

    magfield.XZ.dB2  = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.dB2x = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.dB2y = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    magfield.XZ.dB2z = np.zeros((magfield.z.shape[0],magfield.x.shape[0]))
    
    magfield.YZ.dB2  = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.dB2x = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.dB2y = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))
    magfield.YZ.dB2z = np.zeros((magfield.z.shape[0],magfield.y.shape[0]))

    mu_0 = np.pi*4.e-7
    
    print('First loop')
    for i in range(1,magfield.x.shape[0]-1):
        for j in range(1,magfield.y.shape[0]-1):
            magfield.XY.dB2x[j,i] = \
                 (magfield.XY.Bx[j,i]*(magfield.XY.Bx[j,i+1]-magfield.XY.Bx[j,i-1])/ \
                 (magfield.gstep[2]*magfield.c_omegapi))

            magfield.XY.dB2y[j,i] = \
                (magfield.XY.By[j,i]*(magfield.XY.By[j+1,i]-magfield.XY.By[j-1,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi))

            magfield.XY.dB2z[j,i] = \
                  (magfield.XY.Bz[j,i]*(magfieldsup.XY.Bz[j,i]-magfieldinf.XY.Bz[j,i])/ \
                  (magfield.gstep[0]*magfield.c_omegapi))       
                                                             
            magfield.XY.dB2[j,i] = np.sqrt(magfield.XY.dB2x[j,i]*magfield.XY.dB2x[j,i] \
                                           + magfield.XY.dB2y[j,i]*magfield.XY.dB2y[j,i] \
                                           + magfield.XY.dB2z[j,i]*magfield.XY.dB2z[j,i])  
    print('First loop end - starting second loop')
    print(magfield.XY.dB2[0:20,50])
    print(magfield.c_omegapi)
    print(magfield.gstep)
    
    for i in range(1,magfield.x.shape[0]-1):
        for j in range(1,magfield.z.shape[0]-1):
            magfield.XZ.dB2x[j,i] = \
                (magfield.XZ.Bx[j,i]*(magfield.XZ.Bx[j,i+1]-magfield.XZ.Bx[j,i-1])/ \
                (magfield.gstep[2]*magfield.c_omegapi))

            magfield.XZ.dB2y[j,i] = \
                (magfield.XZ.By[j,i]*(magfieldsup.XZ.By[j,i]-magfieldinf.XZ.By[j,i])/ \
                (magfield.gstep[1]*magfield.c_omegapi))

            magfield.XZ.dB2z[j,i] = \
                (magfield.XZ.Bz[j,i]*(magfield.XZ.Bz[j+1,i]-magfield.XZ.Bz[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi)) 
                                                             
            magfield.XZ.dB2[j,i] = np.sqrt(magfield.XZ.dB2x[j,i]*magfield.XZ.dB2x[j,i] \
                                           + magfield.XZ.dB2y[j,i]*magfield.XZ.dB2y[j,i] \
                                           + magfield.XZ.dB2z[j,i]*magfield.XZ.dB2z[j,i])  
    print('Second loop end - starting third loop')        
    print(magfield.XZ.dB2[0:20,50])

    for i in range(1,magfield.y.shape[0]-1):
        for j in range(1,magfield.z.shape[0]-1):
            magfield.YZ.dB2x[j,i] = \
                (magfield.YZ.Bx[j,i]*(magfieldsup.YZ.Bx[j,i]-magfield.YZ.Bx[j,i])/ \
                (magfield.gstep[2]*magfield.c_omegapi))

            magfield.YZ.dB2y[j,i] = \
                (magfield.YZ.By[j,i]*(magfield.YZ.By[j,i+1]-magfield.YZ.By[j,i-1])/ \
                (magfield.gstep[1]*magfield.c_omegapi))

            magfield.YZ.dB2z[j,i] = \
                (magfield.YZ.Bz[j,i]*(magfield.YZ.Bz[j+1,i]-magfield.YZ.Bz[j-1,i])/ \
                (magfield.gstep[0]*magfield.c_omegapi))      
                                                             
            magfield.YZ.dB2[j,i] = np.sqrt(magfield.YZ.dB2x[j,i]*magfield.YZ.dB2x[j,i] \
                                           + magfield.YZ.dB2y[j,i]*magfield.YZ.dB2y[j,i] \
                                           + magfield.YZ.dB2z[j,i]*magfield.YZ.dB2z[j,i])  
         
    print(magfield.YZ.dB2[0:20,50])  
