# Grid Cell Activity Simulation and Analysis

This repository contains code for simulating and analyzing the correlation of grid cell (across single or multiple module) activity in various configurations [different levels (or absence of) variation in orientation, spacing and phase].

## Files and Directories

- **ipynb/**: Contains Jupyter notebooks for interactive analysis. Core ideas built on ipynb before translating to py.
  - `ActivityDiff.ipynb`: Notebook for analyzing activity differences.
  - `range_analytic_verify.ipynb`: Notebook for verifying range analytics.
- **py/**: Contains Python scripts for simulations and plotting.
  - `Core.py`: Core functions and classes used across the project.
  - `Single_Mod_CorrPlot.py`: Script for plotting correlation of grid cell activity in a single module.
  - `Multi_Mod_CorrPlot.py`: Script for plotting correlation of grid cell activity in multiple modules.
  - `PaperPlot.py`: Script for generating line plots for the paper (Line plots).
  - `PaperPlot2.py`: Second iteration of the plotting script for the paper (Scatter plot).