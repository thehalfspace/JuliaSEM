############################################
### ANALYZE RESULTS FROM THE OUTPUT FILES
############################################

include("output.jl")

include("scripts/earthquake-cycles.jl")
include("scripts/plots.jl")
include("scripts/cumulative-slip.jl")

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

    figname = "/Users/prith/JuliaSEM/plots/test13/stressdrop2.png"
    fig[:savefig](figname, dpi = 300)
    #  figname = string(dir, "/plots", name, "/Vfmax.png")
    #  fig[:savefig](figname, dpi = 300)
end
