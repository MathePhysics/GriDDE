import matplotlib.pyplot as plt
import numpy as np
from Core import *

# Defined in module Core.py, change both here and there for appropriate effect
N = 20
R = 10
num_points = [11, 11, 11]                    # numper of points taken for sampling
# num_points = [21, 21, 21]                    # numper of points taken for sampling

# Globals
arenaSize = 1                       # Size of arena (in meters)
gausswidth = 0.16                    # Width of the gaussian used for each cell

nNeuron = 128                        # To be used in each subplot (length = number of plotlines in each subplot)
plotOrientations = [0, 30]                      # Mean orientation of grid cells for each plot
                                                # length of this = number of plots
nFigs = len(plotOrientations)                   # Number of plots

nModules = [1, 3, 1]                            # nModule for the 3 subplots
gridSpacing = [1, [0.71,1,1.41], 1]             # Mean spacing of grid cells for each module in the subplots
gridOrientation = [0, [0,0,0], 0]               # Mean orientation of grid cells for each module in the subplots
                                                # Nested list for multimodule

spacingStds = [0, 0, 0]                      # Std for sampling grid spacing for the 3 subplots
oriStds = [0, 0, 5]                             # Std for sampling grid orientations for the 3 subplots
nSubplots = len(nModules)                       # Number of subplots

def fig_plotter(nNeuron=nNeuron, plotOrientations=plotOrientations, nModules=nModules, gridSpacing=gridSpacing, spacingStds=spacingStds):
    """
    Plots figures for the paper based on the provided parameters.
    Parameters:
    -----------
    nNeuron : int
        Number of neurons to be used in the plots.
    plotOrientations : list of float
        List of orientations (in degrees) to be used for plotting.
    nModules : list of int
        List of module counts for each subplot.
    gridSpacing : list of float or list of lists of float
        List of grid spacings or list of lists of grid spacings for each module.
    spacingStds : list of float
        List of standard deviations for grid spacings.
    Returns:
    --------
    None
        The function saves the generated plots as PNG and SVG files in the "./Results/" directory.
    Notes:
    ------
    - The function generates subplots for each combination of neuron count and module configuration.
    - Each subplot shows the correlation of activity along a specified direction.
    - The plots are saved with filenames indicating the orientation used.
    """

    Figs, Axs = [],[]
    for _ in range(nFigs):
        fig, axs = plt.subplots(1, 3, figsize=(12,4))
        Figs.append(fig)
        Axs.append(axs)

    for i_subplot in range(nSubplots):
        nModule = nModules[i_subplot]
        spacingStd = spacingStds[i_subplot]
        oriStd = oriStds[i_subplot]

        if nModule == 1:
            orientations =  np.random.normal(gridOrientation[i_subplot], oriStd, nNeuron)
            spacings = np.random.normal(gridSpacing[i_subplot], spacingStd, nNeuron)
            phases = np.random.uniform(-arenaSize, arenaSize, (nNeuron,2))
        else:
            orientations = []
            spacings = []
            phases = []
            for i_module in range(nModule):
                orientations.append( np.random.normal(gridOrientation[i_subplot][i_module], oriStd, nNeuron) )
                spacings.append( np.random.normal(gridSpacing[i_subplot][i_module], spacingStd, nNeuron) )
                phases.append( np.random.uniform(-arenaSize, arenaSize, (nNeuron,2)) )
            orientations = np.concatenate(orientations)
            spacings = np.concatenate(spacings)
            phases = np.concatenate(phases)
        params = list(zip(spacings, phases, orientations))
        Xs, Ys = [], []
        for param in params:
            X,Y = Lattice(*param, N=20)
            Xs.append(X)
            Ys.append(Y)

        def corralong_dir(ax, t=0):
            r = np.linspace(-R,R,num_points[i_subplot])
            theta = np.radians(t)
            x = r*np.cos(theta)
            y = r*np.sin(theta)
            activity_init = GaussLattice(0, 0, Xs, Ys, gausswidth)
            activity = GaussLattice(x, y, Xs, Ys, gausswidth)
            dot_product = np.dot(activity, activity_init)  # Shape (400,)
            norm_activity_init = np.linalg.norm(activity_init)
            norm_activity = np.linalg.norm(activity, axis = 1)
            corr_activity = dot_product/ (norm_activity_init * norm_activity)

            ax.scatter(r, corr_activity, linewidth=2)
            ax.set_title(f"{nModule} modules: " + r"$\Delta \theta$" + f" = {oriStd}" + r"$^{\circ}$, $\Delta \lambda$" + f" = {spacingStd}")
            ax.set_ylim((0,1.1))

            # Set the x and y axis lines to cross at the origin
            ax.spines['left'].set_position('zero')
            ax.spines['bottom'].set_position('zero')

            # Hide the top and right spines (bounding box)
            ax.spines['top'].set_color('none')
            ax.spines['right'].set_color('none')

            # Hide the ticks on the top and right spines
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            # Remove y axis ticks except for 1
            ax.yaxis.set_major_locator(plt.MaxNLocator(1))
            ax.set_xlim(-11,11)
            # Remove the 0 tick in y axis
            yticks = ax.yaxis.get_major_ticks()
            yticks[0].label1.set_visible(False)

            # ax.legend(loc='lower left')

        for i,angle in enumerate(plotOrientations):
            axs = Axs[i]
            fig = Figs[i]
            corralong_dir(axs[i_subplot], angle)


    for i, fig in enumerate(Figs):
        fig.suptitle(f"Orientation = {plotOrientations[i]}", y=0.95)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        fig.savefig(f"./Results/Ori_{plotOrientations[i]}_scatter.png")
        fig.savefig(f"./Results/Ori_{plotOrientations[i]}_scatter.svg")

if __name__ == "__main__":
    fig_plotter()

