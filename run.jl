# Running the simulations from here


include("src/main.jl")

using JLD

s = space_parameters()
tim = time_parameters()
m = medium_properties()
eq = earthquake_parameters()

tic()
FltX, delf5yr, delfsec, Stress, SlipVel, Slip, time_, cca, ccb = 
                                                    @time main(s, tim, m, eq);

println("\n")
elapsed_time = toc()

info("Simulation Complete!")

@save "data/testing/dump01.jld"
