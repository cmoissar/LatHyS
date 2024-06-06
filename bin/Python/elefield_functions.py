import numpy as np
import math as m

from class_data import *
from energy_density import *

def calculate_norm(elefield):

    ## Calculates the magnitude of the electric field
    elefield.E = np.sqrt(np.add(np.add(np.multiply(elefield.Ex,elefield.Ex),
        np.multiply(elefield.Ey,elefield.Ey)),np.multiply(elefield.Ez,elefield.Ez)))
