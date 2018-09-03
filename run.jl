#################################
# Run the simulations from here
#################################

# 1. Go to src/parameters/defaultParameters and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts function


#  using Distributed
#  addprocs(4)

#@everywhere using Distributed
#@everywhere using JLD2
@everywhere using Printf
@everywhere using LinearAlgebra
@everywhere using DelimitedFiles
@everywhere using SharedArrays

@everywhere include("$(@__DIR__)/src/parameters/defaultParameters.jl")	    #	Set Parameters
@everywhere include("$(@__DIR__)/src/setup.jl")

@everywhere P = setParameters(5e3)
@everywhere S = setup(P)

@everywhere include("$(@__DIR__)/src/PCG.jl")               # Preconditioned conjugate gradient to invert matrix
@everywhere include("$(@__DIR__)/src/dtevol.jl")            # compute the next timestep
@everywhere include("$(@__DIR__)/src/NRsearch.jl")          # Newton-rhapson search method to find roots

@everywhere include("$(@__DIR__)/src/main.jl")

simulation_time = @elapsed O = @time main(P, S)

# Save output to file
using Serialization
open("$(@__DIR__)/output/data01.out", "w") do f
    serialize(f,O)
    serialize(f, simulation_time)
end
# Name of the current simulation
#global name = "/dump01"

# Get the current directory for saving figures
#global dir = pwd()


#  simulation_time = @elapsed output = main(parameters(), setup(parameters()))


println("\n")

@info("Simulation Complete!");

# Directory to save the simulation results
#filename = string(dir, "/data", name, ".jld2");

#@save filename output simulation_time name dir 

# Create a new directory to save plots
#mkdir(string(dir, "/plots", name));



