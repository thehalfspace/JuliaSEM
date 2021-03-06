############################################
### ANALYZE RESULTS FROM THE OUTPUT FILES
############################################

include("output.jl")

include("scripts/earthquake-cycles.jl")
include("scripts/plots.jl")
include("scripts/cumulative-slip.jl")

# path to save files
global path = "/Users/prith/JuliaSEM/plots/res04/"

# Deserialize the output
using Serialization
open("data/res04.out") do f
    global O, sim_time, P, S
    O = deserialize(f)
    sim_time = deserialize(f)
    P = deserialize(f)
    S = deserialize(f)
end


#  delfsec, delf5yr = cumSlip(O.Slip, O.SlipVel, O.time_)
delfsec = O.seismic_slip
delf5yr = O.is_slip
delfafter = O.delfafter
stressdrops = O.taubefore-O.tauafter


#  delfafter, stressdrops, tStart, tEnd, vhypo, hypo = Coslip(S, O.Slip, O.SlipVel, O.Stress, O.time_)

Mw, del_sigma = moment_magnitude(P, S, delfafter, stressdrops, O.time_);

function del_sigmaPlot(Mw, del_sigma)

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](Mw, del_sigma, ".")
    ax[:set_xlabel]("Moment Magnitude (Mw)")
    ax[:set_ylabel]("Stress Drops (MPa)")
    ax[:set_title]("Stress Drops vs. Moment Magnitude")
    #  ax[:set_yscale]("log")
    show()

    figname = string(path, "stressdrop.pdf")
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

    figname = string(path, "hypo.pdf")
    fig[:savefig](figname, dpi = 300)
end
