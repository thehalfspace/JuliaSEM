#################################
# Run the simulations from here
#################################

# 1. Go to src/parameters/defaultParameters and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts folder


#using Distributed
#addprocs(4)

#@everywhere using Distributed
@everywhere using JLD2
@everywhere using Printf
@everywhere using LinearAlgebra
@everywhere using DelimitedFiles
@everywhere using SharedArrays

@everywhere include("defaultParameters.jl")	    #	Set Parameters
@everywhere include("setup.jl")

@everywhere P = setParameters(5e3)
@everywhere S = setup(P)

@everywhere include("PCG.jl")               # Preconditioned conjugate gradient to invert matrix
@everywhere include("dtevol.jl")            # compute the next timestep
@everywhere include("NRsearch.jl")          # Newton-rhapson search method to find roots

@everywhere include("main.jl")

simulation_time = @elapsed @time O = main(P, S)

# Name of the current simulation
#global name = "/testparallel"

# Get the current directory for saving figures
#global dir = "" #pwd()


#  simulation_time = @elapsed output = main(parameters(), setup(parameters()))


println("\n")

@info("Simulation Complete!");

# Directory to save the simulation results
#filename = string(dir, "/data", name, ".jld2");

@save "testparallel.jld2" O simulation_time 

# Create a new directory to save plots
#mkdir(string(dir, "/plots", name));
