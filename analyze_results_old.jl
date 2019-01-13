############################################
### ANALYZE RESULTS FROM THE OUTPUT FILES
############################################

include("output_old.jl")


include("scripts_old/scripts/earthquake-cycles.jl")
include("scripts_old/scripts/plots.jl")
include("scripts_old/scripts/cumulative-slip.jl")

# path to save files
global path = "/Users/prith/JuliaSEM/plots/homodc04/"

# Deserialize the output
using Serialization
open("output/wozhi_sims/homodc04.out") do f
    global O, sim_time, P, S
    O = deserialize(f)
    sim_time = deserialize(f)
    P = deserialize(f)
    S = deserialize(f)
end


delfsec, delf5yr = cumSlip(O.Slip, O.SlipVel, O.time_)

delfafter, stressdrops, tStart, tEnd, vhypo, hypo = Coslip(S, O.Slip, O.SlipVel, O.Stress, O.time_)

Mw, del_sigma = moment_magnitude(P, S, O.Slip, O.SlipVel, O.Stress, O.time_);

function del_sigmaPlot(Mw, del_sigma)

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](Mw, del_sigma, ".")
    ax[:set_xlabel]("Moment Magnitude (Mw)")
    ax[:set_ylabel]("Stress Drops (MPa)")
    ax[:set_title]("Stress Drops vs. Moment Magnitude")
    #  ax[:set_yscale]("log")
    show()

    figname = string(path, "stressdrop.png")
    fig[:savefig](figname, dpi = 300)
end

function MwHypoPlot(Mw, hypo)

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](Mw, hypo./1e3, ".")
    ax[:set_xlabel]("Moment Magnitude (Mw)")
    ax[:set_ylabel](" Depth (km)")
    ax[:set_title]("Magnitude vs. Depth")
    show()

    figname = string(path, "hypo.png")
    fig[:savefig](figname, dpi = 300)
end
