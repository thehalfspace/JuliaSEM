This project is out of date, I am maintaining a newer version [here](https://github.com/thehalfspace/eqcycle). 

Spectral element method for earthquake cycle simulations with dynamic treatment of inertial effects.

Adapted from Kaneko et al. (2011), and J.P. Ampuero's [SEMLAB](https://www.mathworks.com/matlabcentral/fileexchange/6154-semlab).

1. change parameters in src/parameters/defaultParameters.jl 
2. edit the initial conditions in src/initialConditions/defaultInitialConditions
3. Change the output file name in run.jl, add number of cpu processors (default 4)
4. run the program in terminal using "julia run.jl"
5. Wait for eternity before the results pop up! (I am still working on making it faster)
