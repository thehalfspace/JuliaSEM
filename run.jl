# Running the simulations from here

include("src/setup.jl")

@time include("src/timesolver.jl")

using JLD
@save "data/testing/test06.jld"
