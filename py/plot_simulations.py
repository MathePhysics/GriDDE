import matplotlib.pyplot as plt
import numpy as np
from Core import Lattice, GaussLattice

ARENA_SIZE = 1 # size of arena (in meters)
R = 1 # maximum distance for correlation 
N_POINTS = 100 # number of points along the line
NEURON_COUNTS = [16, 32, 64] # different numbers of grid cells to simulate
M = 4 # number of modules for case 3

#parameters for each case
case1_params = {"grid_orientation": 0, "ori_std": 10, "grid_spacing": 1, "spacing_std": 0.2}
case2_params = {"grid_orientation": 0, "ori_std": 0, "grid_spacing": 1, "spacing_std": 0}
case3_params = {"orientation_range": (0, 180), "spacing_range": (0.8, 1.2)}

def generate_case3_params(n_neurons: int, orientation_range: tuple, spacing_range: tuple) -> tuple:
    """
    generate grid cell parameters for case 3 (multiple modules with no variability)
    """

    cells_per_module = n_neurons // M
    orientations = []
    spacings = []
    phases = []

    for i in range(M):
        module_orientation = np.random.uniform(*orientation_range)
        module_spacing = np.random.uniform(*spacing_range)
        module_phases = np.random.uniform(-ARENA_SIZE, ARENA_SIZE, (cells_per_module, 2))
        orientations.extend([module_orientation] * cells_per_module)
        spacings.extend([module_spacing] * cells_per_module)
        phases.extend(module_phases)

    # handle remaining cells if n_neurons is not divisible by m
    remaining_cells = n_neurons - M * cells_per_module
    if remaining_cells > 0:
        for _ in range(remaining_cells):
            module_orientation = np.random.uniform(*orientation_range)
            module_spacing = np.random.uniform(*spacing_range)
            module_phase = np.random.uniform(-ARENA_SIZE, ARENA_SIZE, (1, 2))
            orientations.append(module_orientation)
            spacings.append(module_spacing)
            phases.append(module_phase[0])

    return np.array(spacings), np.array(phases), np.array(orientations)

for n_neurons in NEURON_COUNTS:
    print(f"simulating for n_neurons = {n_neurons}")

    # case 1: one module with variability in spacing and orientation
    orientations = np.random.normal(case1_params["grid_orientation"], case1_params["ori_std"], n_neurons)
    spacings = np.random.normal(case1_params["grid_spacing"], case1_params["spacing_std"], n_neurons)
    phases = np.random.uniform(-ARENA_SIZE, ARENA_SIZE, (n_neurons, 2))

    params = list(zip(spacings, phases, orientations))
    X_case1, Y_case1 = [], []
    for param in params:
        X, Y = Lattice(*param, N=10)
        X_case1.append(X)
        Y_case1.append(Y)

    # case 2: one module with no variability
    orientations = np.full(n_neurons, case2_params["grid_orientation"])
    spacings = np.full(n_neurons, case2_params["grid_spacing"])
    phases = np.random.uniform(-ARENA_SIZE, ARENA_SIZE, (n_neurons, 2))

    params = list(zip(spacings, phases, orientations))
    X_case2, Y_case2 = [], []
    for param in params:
        X, Y = Lattice(*param, N=10)
        X_case2.append(X)
        Y_case2.append(Y)

    # case 3: multiple modules with no variability
    spacings, phases, orientations = generate_case3_params(
        n_neurons, case3_params["orientation_range"], case3_params["spacing_range"]
    )

    params = list(zip(spacings, phases, orientations))
    X_case3, Y_case3 = [], []
    for param in params:
        X, Y = Lattice(*param, N=10)
        X_case3.append(X)
        Y_case3.append(Y)

    def compute_and_plot_correlation(case: int, title: str, X: list, Y: list, theta: int):
        """
        compute and plot correlation function for a given theta direction
        """

        r = np.linspace(-R, R, N_POINTS)
        theta_rad = np.radians(theta) 
        x = r * np.cos(theta_rad)
        y = r * np.sin(theta_rad)

        activity_init = GaussLattice(0, 0, X, Y)
        activity = GaussLattice(x, y, X, Y)
        corr_activity = np.sqrt(np.sum(activity * activity_init, axis=1) / n_neurons)

        plt.figure()
        plt.plot(r, corr_activity)
        plt.title(f"{title} (n={n_neurons}, theta={theta}°)")
        plt.xlabel("distance")
        plt.ylabel("correlation")
        plt.savefig(f"Results/plot_cases/correlation_case{case}_n{n_neurons}_theta{theta}.png")
        plt.close()

    for theta in [0, 45, 90]: 
        compute_and_plot_correlation(1, "Case 1: One module with variability", X_case1, Y_case1, theta)
        compute_and_plot_correlation(2, "Case 2: One module with no variability", X_case2, Y_case2, theta)
        compute_and_plot_correlation(3, "Case 3: Multiple modules with no variability", X_case3, Y_case3, theta)
