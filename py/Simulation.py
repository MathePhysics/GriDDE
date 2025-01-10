import matplotlib.pyplot as plt
# from matplotlib.backends.backend_agg import FigureCanvasAgg
# from matplotlib.figure import Figure
# from matplotlib.animation import FuncAnimation  
import numpy as np
import scipy as sp
from tqdm import tqdm
from itertools import product
from Core import *

# Defined in module Core, change both here and there for appropriate effect
N = 50
R = 10
num_points = 400

# Set random seed
# np.random.seed(1)

# Globals
gridOrientation = 0                 # Mean orientation of grid cells (in degrees)
oriStd = 10                         # Std for sampling grid orientation
gridSpacing = 1                     # Mean spacing of grid cells (in meters)    
spacingStd = 0                      # Std for sampling grid spacing
arenaSize = 1                       # Size of arena (in meters)       
gausswidth = 0.3                    # Width of the gaussian used for each cell
nNeurons = 32                       # Sizes of populations of grid cells to be tested
# nSamples = 2                      # Number of independent populations to test per size
# nDecoding = 10                    # Number of random spiking vectors to draw per position
# gridFiring_max_mean = 13          # Mean max firing rate for idealized model
# gridFiring_max_std = 8            # Std for max firing rate
nModules = 4                        # Number of modules 

configstr = f"gridOrientation:{gridOrientation}, oriStd:{oriStd}, gridSpacing:{gridSpacing}, spacingStd:{spacingStd}, arenaSize:{arenaSize}, gausswidth:{gausswidth}, nModules:{nModules} "
meta = { "GRIDCONFIG": configstr }
print(meta)

# Case M1V0 ========== 1 module no variation ====================================
c = 1

# Case M1V1 ========== 1 module with variation ==================================

# Case MnV0 ========== Many module no variation =================================

# Case MnV1 ========== Many module with variation =============================== 

# Randomly generate nNeuron number of grid cells
orientations = np.random.normal(gridOrientation, oriStd, nNeurons)
spacings = np.random.normal(gridSpacing, spacingStd, nNeurons)
phases = np.random.uniform(-arenaSize, arenaSize, (nNeurons,2))
# print(phases)

params = list(zip(spacings, phases, orientations))
Xs, Ys = [], []
for param in params:
    X,Y = Lattice(*param, N=20)
    # plt.scatter(X, Y)
    # plt.title("The distributed lattice points")
    Xs.append(X) 
    Ys.append(Y)
# plt.clf()
# plt.show()
# Xs = np.stack(Xs)
# Ys = np.stack(Ys)
# x = np.array([0,1])
# print(GaussLattice(x,x,Xs,Ys))

def corralong_dir(t=0):
    r = np.linspace(-R,R,num_points)
    theta = np.radians(t)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    # print(x.shape)
    activity_init = GaussLattice(0, 0, Xs, Ys, gausswidth)
    # activity_init2 = np.vstack([activity_init]*num_points)
    activity = GaussLattice(x, y, Xs, Ys, gausswidth)
    # print(np.vstack([activity_init]*num_points).shape)
    # print((activity_init).shape)
    # print((activity).shape)

    # print(activity)
    # print(activity_init)
    dot_product = np.dot(activity, activity_init)  # Shape (400,)
    norm_activity_init = np.linalg.norm(activity_init)
    norm_activity = np.linalg.norm(activity, axis = 1)

    corr_activity = dot_product/ (norm_activity_init * norm_activity)
    # print(dot_product)
    # print(norm_activity)
    # print("hi")
    # print(norm_activity_init)

    # print(corr_activity.shape)
    # return
    # print((activity_init/activity_init).shape)
    # return
    # corr_activity = np.sqrt(activity*activity_init/(nNeurons))
    # corr_activity = np.sqrt(activity@activity_init/(nNeurons**1))
    # corr_activity = np.linalg.norm((activity-activity_init)/nNeurons, axis=1)
    # print(corr_activity.shape)
    # return

    # print((activity@activity_init).shape)

    # fig = Figure(figsize=(5, 4), dpi=100)
    # canvas = FigureCanvasAgg(fig)
    # ax = fig.add_subplot()

    # ax.plot(r,corr_activity )
    # ax.set_title(r"Correlation along $\theta$ = " +f"{t} , n = {nNeurons}")
    # ax.set_xlabel("Distance")
    # ax.set_ylabel("Correlation")
    # ax.set_ylim((0,1.1))
    # fig.savefig(f"Results/plot_cases2/n{nNeurons}_theta{t}_width{gausswidth}.png", metadata=meta)

    plt.plot(r,corr_activity )
    plt.title(r"Correlation along $\theta$ = " +f"{t} , n = {nNeurons}")
    plt.xlabel("Distance")
    plt.ylabel("Correlation")
    plt.ylim((0,1.1))
    plt.savefig(f"Results/plot_cases2/n{nNeurons}_theta{t}_width{gausswidth}.png", metadata=meta)
    plt.close()

angles = list(10*i for i in range(0,18))
for angle in angles:
    corralong_dir(angle)
    # break