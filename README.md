# Bloch modes based calculation of Photonic crystal slab dispersion and Q-factors

This repository provides a basic example of Bloch-mode-based model presented in:
> Ovcharenko, A. I., Blanchard, C., Hugonin, J.-P., & Sauvan, C. (2020). Bound states in the continuum in symmetric and asymmetric photonic crystal slabs. Physical Review B, 101(15), 155303. https://doi.org/10.1103/PhysRevB.101.155303

The model consists of two steps: 
1. Propagative Bloch modes parameters calculation
2. Application of model formulas

The Bloch modes effective indices (essentially, propagation constants, <img src="https://render.githubusercontent.com/render/math?math=\beta">'s, divided by <img src="https://render.githubusercontent.com/render/math?math=k_0">) are defined in an infinite periodic medium, their reflection coefficients -- at the intersect between two media.

![](blochmodes.png)

For instance, in case of two propagative BMs, S-matrix takes the following form

<img src="smatrix.png" width="200">

where 

<img src="https://render.githubusercontent.com/render/math?math=t^'_i"> -- transmission coefficient from the outside space plane wave to ith Bloch mode

<img src="https://render.githubusercontent.com/render/math?math=t_i"> -- ith Bloch mode to plane wave transmission

<img src="https://render.githubusercontent.com/render/math?math=r^'"> -- plane wave to plane wave reflection

<img src="https://render.githubusercontent.com/render/math?math=r_{ii}"> -- ith Bloch mode reflection coefficient

<img src="https://render.githubusercontent.com/render/math?math=r_{ij}"> -- ith to jth Bloch modes cross-reflection coeffi cient

File `RCWA_example.m` provides an example calculation of BMs using RCWA implementation within the `Reticolo` package:
> Hugonin, J. P., & Lalanne, P. (2005). Reticolo software for grating analysis. Institut d'Optique Graduate School. Retrieved from https://www.lp2n.institutoptique.fr/equipes-de-recherche-du-lp2n/light-complex-nanostructures

Use this link to download the package, as well as its documentation.

File `BICModel_example.m` uses `BICModel.m` function to calculate dispersion and Q-factor of a PhC mode.
