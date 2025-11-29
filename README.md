# Grid cell Decorrelation for Distance Encoding (GRIDDE)

This repository contains code for simulating and analyzing the (De-)correlation of grid cell activity based on the following article [(arXiv)](https://arxiv.org/abs/2511.08292):\
Distance by de-correlation: Computing distance with heterogeneous grid cells \
Pritipriya Dasbehera, Akshunna S. Dogra and William T. Redman

## Files and Directories

- **ipynb/**: Contains Jupyter notebooks for interactive analysis. Core ideas built on ipynb before translating to py.
  - `ActivityDiff.ipynb`: Notebook for analyzing activity differences.
  - `range_analytic_verify.ipynb`: Notebook for verifying range analytics.
- **py/**: Contains Python scripts for simulations and plotting.
  - `Core.py`: Core functions and classes used across the project.
  - `Single_Mod_CorrPlot.py`: Script for plotting correlation of grid cell activity in a single module.
  - `Multi_Mod_CorrPlot.py`: Script for plotting correlation of grid cell activity in multiple modules.
  - `PaperPlot.py`: Script for generating line plots. (v2 generalises to many neuron counts)
  - `PaperPlot2.py`: Second iteration of the plotting script (Scatter plot). (v2 generalises to many neuron counts)
- **theoretical_plots/powerspectrum.jl**: Contains Julia code to verify the theoretical arguments for 1D case presented in the paper.
