import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation  
import numpy as np
import scipy as sp
from tqdm import tqdm
from itertools import product
from Core import *

# Set random seed
# np.random.seed(1)

# Globals
gridOrientation = 0                 # Mean orientation of grid cells (in degrees)
oriStd = 10                         # Std for sampling grid orientation
gridSpacing = 1                     # Mean spacing of grid cells (in meters)    
spacingStd = 0.2                    # Std for sampling grid spacing
arenaSize = 1                       # Size of arena (in meters)       
nNeurons = 2                        # Sizes of populations of grid cells to be tested
# nSamples = 2                      # Number of independent populations to test per size
# nDecoding = 10                    # Number of random spiking vectors to draw per position
# gridFiring_max_mean = 13          # Mean max firing rate for idealized model
# gridFiring_max_std = 8            # Std for max firing rate

# Randomly generate nNeuron number of grid cells
orientations = np.random.normal(gridOrientation, oriStd, nNeurons)
spacings = np.random.normal(gridSpacing, spacingStd, nNeurons)
phases = np.random.uniform(-arenaSize, arenaSize, (nNeurons,2))
print(phases)

params = list(zip(spacings, phases, orientations))
Xs, Ys = [], []
for param in params:
    X,Y = Lattice(*param, N=10)
    plt.scatter(X, Y)
    plt.title("The distributed lattice points")
    Xs.append(X) 
    Ys.append(Y)
plt.show()
# Xs = np.stack(Xs)
# Ys = np.stack(Ys)
# x = np.array([0,1])
# print(GaussLattice(x,x,Xs,Ys))

def corralong_dir(t=0):
    r = np.linspace(-R,R,num_points)
    theta = np.radians(t)
    x = r*np.cos(theta)
    y = r*np.sin(theta)

    activity_init = GaussLattice(0, 0, Xs, Ys)
    activity = GaussLattice(x, y, Xs, Ys)
    corr_activity = np.sqrt(activity@activity_init/nNeurons)
    # print((activity@activity_init).shape)
    plt.plot(r,corr_activity )
    plt.title(f"Correlation along theta = {t} deg direction")
    plt.xlabel("Distance")
    plt.ylabel("Correlation")
    plt.show()

corralong_dir(45)