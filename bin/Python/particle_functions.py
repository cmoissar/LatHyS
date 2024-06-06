import numpy as np
import math as m
from class_data import *
from energy_density import *

def calculate_velocity(particle):

    particle.V = np.sqrt(np.add(np.add(np.multiply(particle.Vx,particle.Vx),
        np.multiply(particle.Vy,particle.Vy)),np.multiply(particle.Vz,particle.Vz)))

def calc_thermal_pressure(data_particle,calc_grad=0) :
    # Returns the thermal pressure (in nPa)
    print('Calculate Pth')
    k_B = 1.3806488e-23
    # T is given in eV
    if calc_grad == 1 :
        data_particle.Pth = np.multiply(data_particle.n,data_particle.T)*k_B*1.e15*11605 
    else :
        data_particle.Pth = np.multiply(data_particle.n,data_particle.T)*k_B*1.e15*11605                               


def calculate_thermal_pressure_gradient(particle,particlesup,particleinf):

    calc_thermal_pressure(particle.XY,calc_grad=1)
    calc_thermal_pressure(particle.XZ,calc_grad=1)
    calc_thermal_pressure(particle.YZ,calc_grad=1)
    
    calc_thermal_pressure(particlesup.XY,calc_grad=1)
    calc_thermal_pressure(particlesup.XZ,calc_grad=1)
    calc_thermal_pressure(particlesup.YZ,calc_grad=1)    

    calc_thermal_pressure(particleinf.XY,calc_grad=1)
    calc_thermal_pressure(particleinf.XZ,calc_grad=1)
    calc_thermal_pressure(particleinf.YZ,calc_grad=1)

    particle.XY.grad_Pth   = np.zeros((particle.y.shape[0],particle.x.shape[0]))
    particle.XY.grad_Pthx  = np.zeros((particle.y.shape[0],particle.x.shape[0]))
    particle.XY.grad_Pthy  = np.zeros((particle.y.shape[0],particle.x.shape[0]))
    particle.XY.grad_Pthz  = np.zeros((particle.y.shape[0],particle.x.shape[0]))
    
    particle.XZ.grad_Pth   = np.zeros((particle.z.shape[0],particle.x.shape[0]))
    particle.XZ.grad_Pthx  = np.zeros((particle.z.shape[0],particle.x.shape[0]))
    particle.XZ.grad_Pthy  = np.zeros((particle.z.shape[0],particle.x.shape[0]))
    particle.XZ.grad_Pthz  = np.zeros((particle.z.shape[0],particle.x.shape[0]))
        
    particle.YZ.grad_Pth   = np.zeros((particle.z.shape[0],particle.y.shape[0]))
    particle.YZ.grad_Pthx  = np.zeros((particle.z.shape[0],particle.y.shape[0]))
    particle.YZ.grad_Pthy  = np.zeros((particle.z.shape[0],particle.y.shape[0]))
    particle.YZ.grad_Pthz  = np.zeros((particle.z.shape[0],particle.y.shape[0]))

    for i in range(1,particle.x.shape[0]-1):
        for j in range(1,particle.y.shape[0]-1):

            particle.XY.grad_Pthx[j,i] = (particle.XY.Pth[j,i-1] - particle.XY.Pth[j,i+1])/ \
                (particle.gstep[2]*particle.c_omegapi*2.)
            particle.XY.grad_Pthy[j,i] = (particle.XY.Pth[j-1,i] - particle.XY.Pth[j+1,i])/ \
                (particle.gstep[1]*particle.c_omegapi*2.)
            particle.XY.grad_Pthz[j,i] = (particlesup.XY.Pth[j,i] - particleinf.XY.Pth[j,i])/ \
                (particle.gstep[0]*particle.c_omegapi*2.)

            particle.XY.grad_Pth[j,i] = np.sqrt(particle.XY.grad_Pthx[j,i]*particle.XY.grad_Pthx[j,i] +
                                                 particle.XY.grad_Pthy[j,i]*particle.XY.grad_Pthy[j,i] +
                                                 particle.XY.grad_Pthz[j,i]*particle.XY.grad_Pthz[j,i])

    for i in range(1,particle.x.shape[0]-1):
        for j in range(1,particle.z.shape[0]-1):

            particle.XZ.grad_Pthx[j,i] = (particle.XZ.Pth[j,i-1] - particle.XZ.Pth[j,i+1])/ \
                (particle.gstep[2]*particle.c_omegapi*2.)
            particle.XZ.grad_Pthy[j,i] = (particlesup.XZ.Pth[j,i] - particleinf.XZ.Pth[j,i])/ \
                (particle.gstep[1]*particle.c_omegapi*2.)
            particle.XZ.grad_Pthz[j,i] = (particle.XZ.Pth[j+1,i] - particle.XZ.Pth[j-1,i])/ \
                (particle.gstep[0]*particle.c_omegapi*2.)

            particle.XZ.grad_Pth[j,i] = np.sqrt(particle.XZ.grad_Pthx[j,i]*particle.XZ.grad_Pthx[j,i] +
                                                 particle.XZ.grad_Pthy[j,i]*particle.XZ.grad_Pthy[j,i] +
                                                 particle.XZ.grad_Pthz[j,i]*particle.XZ.grad_Pthz[j,i])

    for i in range(1,particle.y.shape[0]-1):
        for j in range(1,particle.z.shape[0]-1):

            particle.YZ.grad_Pthx[j,i] = (particlesup.YZ.Pth[j,i] - particleinf.YZ.Pth[j,i])/ \
                (particle.gstep[2]*particle.c_omegapi*2.)
            particle.YZ.grad_Pthy[j,i] = (particle.YZ.Pth[j,i-1] - particle.YZ.Pth[j,i+1])/ \
                (particle.gstep[1]*particle.c_omegapi*2.)
            particle.YZ.grad_Pthz[j,i] = (particle.YZ.Pth[j+1,i] - particle.YZ.Pth[j-1,i])/ \
                (particle.gstep[0]*particle.c_omegapi*2.)

            particle.YZ.grad_Pth[j,i] = np.sqrt(particle.YZ.grad_Pthx[j,i]*particle.YZ.grad_Pthx[j,i] +
                                                 particle.YZ.grad_Pthy[j,i]*particle.YZ.grad_Pthy[j,i] +
                                                 particle.YZ.grad_Pthz[j,i]*particle.YZ.grad_Pthz[j,i])

