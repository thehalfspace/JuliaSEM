#################################
# READ OUTPUT FROM SIMULATION
#################################

mutable struct results
    Stress::Array{Float64,2}
    SlipVel::Array{Float64,2}
    Slip::Array{Float64,2}
    time_::Array{Float64}
end

using Serialization
open("output/data02.out") do f
    global Op, sim_time
    Op = deserialize(f)
    sim_time = deserialize(f)
end


using Distributed

include("src/parameters/defaultParameters.jl")
include("src/setup.jl")

P = setParameters(6e3)
S = setup(P)



# Cumulative Slip Plot
include("scripts/cumulative-slip.jl")
include("scripts/plots.jl")

delfsec, delf5yr = cumSlip(Op.Slip, Op.SlipVel, Op.time_)

cumSlipPlot(delfsec, delf5yr, S.FltX)
