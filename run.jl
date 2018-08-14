#################################
# Run the simulations from here
#################################

# 1. Go to src/parameters/defaultParameters and change as needed
# 2. Go to src/initialConditions/defaultInitialConditions and change as needed
# 3. Change the name of the simulation in this file
# 4. Run the simulation from terminal. (julia run.jl)
# 5. Plot results from the scripts function

# Name of the current simulation
global name = "/dump02"

# Get the current directory for saving figures
global dir = pwd()

# include the main function
include(string(dir, "/src/main.jl"))

using JLD2

s = space_parameters()
tim = time_parameters()
m = medium_properties()
eq = earthquake_parameters()

elapsed_time = @elapsed FltX, delf5yr, delfsec, Stress, SlipVel, Slip, time_, cca, ccb = 
                                                    @time main(s, tim, m, eq);

println("\n")

@info("Simulation Complete!")

# Directory to save the simulation results
#filename = string(dir, "/data", name, ".jld")

#@save filename

# Create a new directory to save plots
#mkdir(string(dir, "/plots/", name))
