"""
This module is designed to simulate and analyze the correlation of grid cell activity in a SINGLE MODULE.
written by @Pritipriya_dasbehera
"""

import matplotlib.pyplot as plt
import numpy as np
from Core import *

# Defined in module Core.py, change both here and there for appropriate effect
N = 10
R = 3
num_points = 100

# Globals
gridOrientation = 0                 # Mean orientation of grid cells (in degrees)
oriStd = 10                         # Std for sampling grid orientation
gridSpacing = 1                     # Mean spacing of grid cells (in meters)    
spacingStd = 0                      # Std for sampling grid spacing
arenaSize = 1                       # Size of arena (in meters)       
gausswidth = 0.3                    # Width of the gaussian used for each cell
nNeurons = 32                       # Sizes of populations of grid cells to be tested

configstr = f"gridOrientation:{gridOrientation}, oriStd:{oriStd}, gridSpacing:{gridSpacing}, spacingStd:{spacingStd}, arenaSize:{arenaSize}, gausswidth:{gausswidth}"
meta = { "GRIDCONFIG": configstr }
print(meta)

def multicorr_fig(axs: list[plt.Axes], oriStd = oriStd, spacingStd=spacingStd, gausswidth=gausswidth, nNeurons=nNeurons):
    ''' Multicell correlation plot generator for single module'''

    orientations = np.random.normal(gridOrientation, oriStd, nNeurons)
    spacings = np.random.normal(gridSpacing, spacingStd, nNeurons)
    phases = np.random.uniform(-arenaSize, arenaSize, (nNeurons,2))

    params = list(zip(spacings, phases, orientations))
    Xs, Ys = [], []
    for param in params:
        X,Y = Lattice(*param, N=20)
        Xs.append(X) 
        Ys.append(Y)

    # def plot_sample_ratefield():
    #     r = np.linspace(-R,R,200)
    #     x,y = np.meshgrid(r, r)
    #     activity = GaussLattice(x,y, Xs, Ys, gausswidth)
    #     plt.imshow(activity)
    #     plt.show()
    #     plt.close()
    #     print(np.min(activity))

    # print("plotting now")
    # plot_sample_ratefield()
    
    def corralong_dir(ax, t=0):
        r = np.linspace(-R,R,num_points)
        theta = np.radians(t)
        x = r*np.cos(theta)
        y = r*np.sin(theta)
        activity_init = GaussLattice(0, 0, Xs, Ys, gausswidth)
        activity = GaussLattice(x, y, Xs, Ys, gausswidth)
        dot_product = np.dot(activity, activity_init)  
        norm_activity_init = np.linalg.norm(activity_init)
        norm_activity = np.linalg.norm(activity, axis = 1)
        corr_activity = dot_product/ (norm_activity_init * norm_activity)
        
        ax.plot(r,corr_activity )
        # ax.set_title(r"$\theta$ = " +f"{t} , n = {nNeurons}, \n std={spacingStd}, width={gausswidth}", wrap=True)
        ax.set_title(f"std={spacingStd}, width={gausswidth}", wrap=True)
        # ax.set_xlabel("Distance")
        # ax.set_ylabel("Correlation")
        ax.set_ylim((0,1.1))

    global angles
    # angles = list(10*i for i in range(0,18))
    angles = np.linspace(0,180,len(axs))
    for i,angle in enumerate(angles):
        corralong_dir(axs[i], angle)

def multicorr_plotter():
    s = 4
    fig, axs = plt.subplots(nrows=9,ncols=3, figsize=(2*s,3.5*s))

    spacingStds = [0,0.1,0.2]

    for i in range(3):
        multicorr_fig( axs[i], spacingStd=spacingStds[i], gausswidth=0.2 )
    for i in range(3):
        multicorr_fig( axs[3+i], spacingStd=spacingStds[i], gausswidth=0.25 )
    for i in range(3):
        multicorr_fig( axs[6+i], spacingStd=spacingStds[i], gausswidth=0.3 )
    
    # line = plt.Line2D([0, 1], [0.35, 0.35], color="black", linewidth=1, transform=fig.transFigure, linestyle="--")
    # fig.add_artist(line)

    fig.suptitle(f"nNeurons = {nNeurons} OriStd = {oriStd}\n \n", fontsize='xx-large')
    plt.tight_layout(rect=[0, 0, 0.95, 1])

    pos1 = axs[2, 2].get_position().y0
    fig.text(0.95, pos1+0.1 , "width = 0.2", rotation=90, fontsize=14)
    fig.add_artist(plt.Line2D([0, 1], [pos1-0.017, pos1-0.017], linewidth=1, linestyle='--', color='black'))

    pos2 = axs[5, 2].get_position().y0
    fig.text(0.95, pos2+0.1 , "width = 0.3", rotation=90, fontsize=14)
    fig.add_artist(plt.Line2D([0, 1], [pos2-0.017, pos2-0.017], linewidth=1, linestyle='--', color='black'))

    pos3 = axs[8, 2].get_position().y0
    fig.text(0.95, pos3+0.1 , "width = 0.4", rotation=90, fontsize=14)
    # fig.add_artist(plt.Line2D([0, 1], [pos3-0.01, pos3-0.01], linewidth=1, linestyle='--', color='black'))

    # fig.text(0.95, 0.5, "width = 0.2", ha='center', va='center', rotation=90, fontsize=14)
    # fig.text(0.95, 0.5, "width = 0.2", ha='center', va='center', rotation=90, fontsize=14)


    
    for i in range(3):
        pos = axs[0, i].get_position()  
        centerpos = (pos.x0 + pos.x1)/2
        print(centerpos)
        fig.text( centerpos, pos.y1 + 0.03, r"$\theta$ = " +f"{angles[i]}", ha='center', va='center', fontsize=14, fontweight='bold')

    plt.savefig(f"Results/multicell_correlation_n_{nNeurons}_ori_{oriStd}.png")
    # plt.show()


multicorr_plotter()