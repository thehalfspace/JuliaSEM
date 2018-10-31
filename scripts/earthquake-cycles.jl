#################################
# MODULE FOR SOME CALCULATIONS
# FROM SIMULATION OUTPUT
#################################

using StatsBase
using PyPlot

PyPlot.matplotlib[:rc]("patch.force_edgecolor=true")

#-------------------------------
# Compute hypocenter locations
#------------------------------
function plotHypo(S, Slip, SlipVel, Stress, time_)

    delfafter, stressdrops, tStart, tEnd, vhypo, hypo = 
                                        Coslip(S, Slip, SlipVel, Stress, time_)

    # Plot hypocenter
    hist = fit(Histogram, hypo./1e3, nbins = 10)

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:barh](hist.edges[1][1:end-1], hist.weights, 0.5)
    ax[:plot](collect(1:40), -8*ones(40), "--", label="Fault Zone Depth")
    ax[:set_xlabel]("Number of Earthquakes")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Magnitude-frequency distribution")
    ax[:legend](loc="upper right")
    show()

    figname = "/Users/prith/JuliaSEM/plots/test11/hypocenter.png"
    fig[:savefig](figname, dpi = 300)
    #  figname = string(dir, "/plots", name, "/mfd.png")
    #  ax[:savefig](figname, dpi = 300)

end

#.................................................
# Compute the final Coseismic slip for each event
#.................................................
function Coslip(S, Slip, SlipVel, Stress, time_=zeros(1000000))

    Vfmax = maximum(SlipVel, dims = 1)[:]

    delfafter::Array{Float64,2} = zeros(size(Slip))
    tStart::Array{Float64} = zeros(size(Slip[1,:]))
    tEnd::Array{Float64} = zeros(size(Slip[1,:]))

    taubefore::Array{Float64,2} = zeros(size(Slip))
    tauafter::Array{Float64,2} = zeros(size(Slip))
    
    hypo::Array{Float64} =  zeros(size(Slip[1,:]))   # Hypocenter
    vhypo::Array{Float64} = zeros(size(Slip[1,:]))   # Velocity at hypocenter

    Vthres = 0.001 # event threshold
    slipstart = 0
    it = 1; it2 = 1
    delfref = zeros(size(Slip[:,1]))

    for i = 1:length(Slip[1,:])

        # Start of each event
        if Vfmax[i] > 1.01*Vthres && slipstart == 0
            delfref = Slip[:,i]
            slipstart = 1
            tStart[it2] = time_[i]
            
            taubefore[:,it2] = Stress[:,i]
            vhypo[it2], indx = findmax(SlipVel[:,i])

            hypo[it2] = S.FltX[indx]

            it2 = it2+1
        end

        # End of each event
        if Vfmax[i] < 0.99*Vthres && slipstart == 1
            delfafter[:,it] = Slip[:,i] - delfref
            tauafter[:,it] = Stress[:,i]
            tEnd[it] = time_[i]
            slipstart = 0
            it = it + 1
        end
    end

    return delfafter[:,1:it-1], (taubefore-tauafter)[:,1:it-1], tStart[1:it2-1], tEnd[1:it-1], vhypo[1:it2-1], hypo[1:it2-1]
end

#..........................................................
# Compute the moment magnitude:
#       Assumed the rupture area to be square; the rupture
#       dimension along depth is the same as the rupture
#       dimension perpendicular to the plane
#..........................................................
function moment_magnitude(P, S, Slip, SlipVel, Stress, time_)

    # Final coseismic slip of each earthquake
    delfafter, stressdrops = Coslip(S, Slip, SlipVel, Stress, time_)

    iter = length(delfafter[1,:])
    mu = P.rho1*P.vs1^2
    seismic_moment = zeros(iter)
    temp_sigma = 0
    iter2 = 1 

    del_sigma = zeros(iter)
    
    dx = diff(S.FltX)

    for i = 1:iter
        
        # slip threshold = 1% of maximum slip
        slip_thres = 0.10*maximum(delfafter[:,i])

        # area = slip*(rupture dimension along depth)
        # zdim = rupture along z dimension = depth rupture dimension
        area = 0; zdim = 0; temp_sigma = 0

        for j = 1:P.FltNglob
            if delfafter[j,i] >= slip_thres
                area = area + delfafter[j,i]*dx[j-1]
                zdim = zdim + dx[j-1]

                # Avg. stress drops along rupture area
                temp_sigma = temp_sigma + stressdrops[j,i]*dx[j-1]
            end
        end
        
        seismic_moment[i] = mu*area*zdim
        del_sigma[i] = temp_sigma


    end
    seismic_moment = filter!(x->x!=0, seismic_moment)
    del_sigma = filter!(x->x!=0, del_sigma)
    Mw = (2/3)*log10.(seismic_moment.*1e7) .- 10.7

    return Mw, del_sigma
end


#...........
# Plot MFD
#...........
function MwPlot(Mw)

    hist = fit(Histogram, Mw, nbins = 20)

    # Cumulative
    cum = cumsum(hist.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    #  ax[:plot](hist.edges[1][1:end-1], hist.weights, ".", label="Non-cumulative")
    ax[:plot](hist.edges[1][1:end-1], cum, ".", label="Cumulative")
    ax[:set_xlabel]("Moment Magnitude (Mw)")
    ax[:set_ylabel]("Number of Earthquakes")
    ax[:set_yscale]("log")
    ax[:set_title]("Magnitude-frequency distribution")
    ax[:legend](loc="upper right")
    show()

    figname = "/Users/prith/JuliaSEM/plots/test11/mfd.png"
    fig[:savefig](figname, dpi = 300)
    #  figname = string(dir, "/plots", name, "/mfd.png")
    #  ax[:savefig](figname, dpi = 300)
end


#.................................
# Plot earthquake catalog
# (Earthquake magnitude with time)
#.................................
function eq_catalog(Mw, t_catalog, yr2sec)

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:scatter](t_catalog./yr2sec, Mw, s= 30, marker=".")
    ax[:set_xlabel]("Time (yrs)")
    ax[:set_ylabel]("Moment Magnitude (Mw)")
    ax[:set_title]("Earthquake Catalogue")
    show()

    figname = "/Users/prith/JuliaSEM/plots/test11/catalogue.png"
    fig[:savefig](figname, dpi = 300)
    #  figname = string(dir, "/plots", name, "catalogue.png")
    #  fig[:savefig](figname, dpi = 300)
end
