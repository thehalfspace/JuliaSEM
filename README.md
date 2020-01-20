This project is out of data, I am maintaining a newer version [here](https://github.com/thehalfspace/eqcycle). 

Spectral element method for earthquake cycle simulations with dynamic treatment of inertial effects.

Written in Julia 1.0. Syntax is very similar to Matlab but it is much faster.

Adapted from Kaneko et al. (2011), and J.P. Ampuero's [SEMLAB](https://www.mathworks.com/matlabcentral/fileexchange/6154-semlab).

This project is incomplete, work in progress. To run the program 

1. change parameters in src/parameters/defaultParameters.jl 
2. edit the initial conditions in src/initialConditions/defaultInitialConditions
3. Change the output file name in run.jl, add number of cpu processors (default 4)
4. run the program in terminal using "julia run.jl"
5. Wait for eternity before the results pop up! (I am still working on making it faster)
