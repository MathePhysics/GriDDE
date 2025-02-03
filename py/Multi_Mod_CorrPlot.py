"""
This module is designed to simulate and analyze the correlation of grid cell activity in a MULTIPLE MODULES.
written by @Pritipriya_dasbehera
"""

import matplotlib.pyplot as plt
import numpy as np
from Core import *

# Defined in module Core.py, change both here and there for appropriate effect
N = 10
R = 3
num_points = 100

# Globals (define a list for data at each module)
gridOrientation = [0]       # Mean orientation of grid cells for each module (in degrees)
oriStd = 0                          # Std for sampling grid orientation (used in multicorr_fig)
oriStds = [0, 1, 2]                # Std for sampling grid orientations (used in multicorr_plotter)

gridSpacing = [1]       # Mean spacing of grid cells (in meters)    
spacingStd = 0                      # Std for sampling grid spacing (used in multicorr_fig)
spacingStds = [0, 0.03, 0.05]         # Std for sampling grid spacing (used in multicorr_plotter)

arenaSize = 1                       # Size of arena (in meters)       
gausswidth = 0.2                    # Width of the gaussian used for each cell
nNeurons = [32, 64, 256]                       # Sizes of populations of grid cells in each module to be tested 
                                    # [list] to re-run everything and experiment at different n values
nModules = len(gridOrientation)     # Number of modules 

configstr = f"gridOrientation:{gridOrientation}, oriStd:{oriStd}, gridSpacing:{gridSpacing}, spacingStd:{spacingStd}, arenaSize:{arenaSize}, gausswidth:{gausswidth}, nModules:{nModules} "
meta = { "GRIDCONFIG": configstr }
print(meta)

def multicorr_fig(axs: list[plt.Axes], oriStd = oriStd, spacingStd=spacingStd, gausswidth=gausswidth, nNeurons=nNeurons):
    ''' Multicell correlation plot generator for multiple(or single) modules'''

    orientations = []
    spacings = []
    phases = []
    for i in range(nModules):
        orientations.append( np.random.normal(gridOrientation[i], oriStd, nNeurons) )
        spacings.append( np.random.normal(gridSpacing[i], spacingStd, nNeurons) )
        phases.append( np.random.uniform(-arenaSize, arenaSize, (nNeurons,2)) )

    orientations = np.concatenate(orientations)
    spacings = np.concatenate(spacings)
    phases = np.concatenate(phases)

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
        dot_product = np.dot(activity, activity_init)  # Shape (400,)
        norm_activity_init = np.linalg.norm(activity_init)
        norm_activity = np.linalg.norm(activity, axis = 1)
        corr_activity = dot_product/ (norm_activity_init * norm_activity)
        
        ax.plot(r,corr_activity )
        # ax.set_title(r"$\theta$ = " +f"{t} , n = {nNeurons}, \n std={spacingStd}, width={gausswidth}", wrap=True)
        ax.set_title(f"std = {spacingStd}", wrap=True)
        # ax.set_xlabel("Distance")
        # ax.set_ylabel("Correlation")
        ax.set_ylim((0,1.1))

    global angles
    # angles = list(10*i for i in range(0,18))
    angles = np.linspace(0,30,len(axs))
    for i,angle in enumerate(angles):
        corralong_dir(axs[i], angle)

def multicorr_plotter(nNeurons, spacingStds, oriStds):
    s = 4
    fig, axs = plt.subplots(nrows=9,ncols=3, figsize=(2*s,3.5*s))

    for i in range(3):
        multicorr_fig( axs[i], spacingStd=spacingStds[i], oriStd=oriStds[i], nNeurons=nNeurons )
    for i in range(3):
        multicorr_fig( axs[3+i], spacingStd=spacingStds[i], oriStd=oriStds[i], nNeurons=nNeurons )
    for i in range(3):
        multicorr_fig( axs[6+i], spacingStd=spacingStds[i], oriStd=oriStds[i], nNeurons=nNeurons )
    
    # line = plt.Line2D([0, 1], [0.35, 0.35], color="black", linewidth=1, transform=fig.transFigure, linestyle="--")
    # fig.add_artist(line)

    fig.suptitle(f"nNeurons X nModules = {nNeurons} X {nModules} \n orientations = {gridOrientation} , spacings = {gridSpacing} \n \n", fontsize='xx-large')
    plt.tight_layout(rect=[0, 0, 0.95, 1])

    posdisplacement = 0.02
    pos1 = axs[2, 2].get_position().y0
    fig.text(0.95, pos1+0.1 , f"oriStd = {oriStds[0]}", rotation=90, fontsize=14)
    fig.add_artist(plt.Line2D([0, 1], [pos1-posdisplacement, pos1-posdisplacement], linewidth=1, linestyle='--', color='black'))

    pos2 = axs[5, 2].get_position().y0
    fig.text(0.95, pos2+0.1 , f"oriStd = {oriStds[1]}", rotation=90, fontsize=14)
    fig.add_artist(plt.Line2D([0, 1], [pos2-posdisplacement, pos2-posdisplacement], linewidth=1, linestyle='--', color='black'))

    pos3 = axs[8, 2].get_position().y0
    fig.text(0.95, pos3+0.1 , f"oriStd = {oriStds[2]}", rotation=90, fontsize=14)
    # fig.add_artist(plt.Line2D([0, 1], [pos3-0.01, pos3-0.01], linewidth=1, linestyle='--', color='black'))

    # fig.text(0.95, 0.5, "width = 0.2", ha='center', va='center', rotation=90, fontsize=14)
    # fig.text(0.95, 0.5, "width = 0.2", ha='center', va='center', rotation=90, fontsize=14)


    
    for i in range(3):
        pos = axs[0, i].get_position()  
        centerpos = (pos.x0 + pos.x1)/2
        fig.text( centerpos, pos.y1 + 0.03, r"$\theta$ = " +f"{angles[i]}", ha='center', va='center', fontsize=14, fontweight='bold')

    # plt.savefig(f"Results/multimod_correlation_root2spacing_n_{nNeurons}_m_{nModules}.png")
    plt.savefig(f"Results/singlemod_correlation_n_{nNeurons}.png")
    # plt.show()

if isinstance(nNeurons, int):
    multicorr_plotter(nNeurons)
else:
    for n in nNeurons:
        multicorr_plotter(n,spacingStds, oriStds)
        print(f"nNeuron = {n} Done")