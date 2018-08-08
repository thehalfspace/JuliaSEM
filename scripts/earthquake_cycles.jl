#################################
# MODULE FOR SOME CALCULATIONS
# FROM SIMULATION OUTPUT
#################################

using StatsBase
using PyPlot

#............................................
# Compute maximum slip rate at each timestep
#............................................
function Vfmax(SlipVel)
    return maximum(SlipVel, 1)
end

#.................................................
# Compute the final Coseismic slip for each event
#.................................................
function Coslip(Slip, SlipVel)

    Vfmax = maximum(SlipVel, 1)

    delfafter::Array{Float64,2} = zeros(Slip)
    Vthres = 0.01 # event threshold
    slipstart = 0
    it = 1
    delfref = zeros(Slip[:,1])

    for i = 1:length(Slip[1,:])

        # Start of each event
        if Vfmax[i] > 1.01*Vthres && slipstart == 0
            delfref = Slip[:,i]
            slipstart = 1
        end

        # End of each event
        if Vfmax[i] < 0.99*Vthres && slipstart == 1
            delfafter[:,it] = Slip[:,i] - delfref
            slipstart = 0
            it = it + 1
        end
    end

    return delfafter[:,1:it-1]
end

#..........................................................
# Compute the moment magnitude:
#       Assumed the rupture area to be square; the rupture
#       dimension along depth is the same as the rupture
#       dimension perpendicular to the plane
#..........................................................
function moment_magnitude(s, m, Slip, SlipVel)

    # Final coseismic slip of each earthquake
    delfafter = Coslip(Slip, SlipVel)

    iter = length(delfafter[1,:])
    
    moment = zeros(iter)

    for i = 1:iter
        
        # slip threshold = 1% of maximum slip
        slip_thres = 0.01*maximum(delfafter[:,i])

        # area = slip*(rupture dimension along depth)
        # zdim = rupture along z dimension = depth rupture dimension
        area = 0; zdim = 0

        for j = 1:s.FltNglob
            if delfafter[j,i] >= slip_thres
                area = area + delfafter[j,i]*s.dxe
                zdim = zdim + s.dxe
            end
        end

        moment[i] = m.mu[1,1]*area*zdim

    end

    Mw = (2/3)*log10.(moment.*1e7) - 10.7

    return Mw
end


#...........
# Plot MFD
#...........
function MwPlot(Mw)

    hist = fit(Histogram, Mw, nbins = 10)

    # Cumulative
    cum = cumsum(hist.weights[end:-1:1])[end:-1:1]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](hist.edges[1][1:end-1], log10.(hist.weights), ".", label="Non-cumulative")
    ax[:plot](hist.edges[1][1:end-1], log10.(cum), ".", label="Cumulative")
    ax[:set_xlabel]("Moment Magnitude (Mw)")
    ax[:set_ylabel]("Number of Earthquakes (log10)")
    ax[:set_title]("Magnitude-frequency distribution")
    ax[:legend](loc="upper right")
    show()
end
