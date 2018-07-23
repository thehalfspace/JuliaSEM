###################################
# MAGNITUDE-FREQUENCY DISTRIBUTION
# FOR EARTHQUAKE CYCLE SIMULATIONS
###################################

# Calculate magnitude
function moment_magnitude(G, Coslip, dummy_it, SlipVel, dxe)

    Mw = zeros(length(dummy_it)-1)

    for ev = 2:length(dummy_it)
        
        iter = Int(dummy_it[ev])

        # find where slip velocity is > 1e-3
        arr = find(SlipVel[:, iter] .>= 1e-3)

        # Integrate slip with rupture area to get moment
        Mo = G*sum(abs.(Coslip[arr,ev]).*dxe^2)

        # Moment magnitude (1 N-m = 1e7 dyne-cm)
        Mw[ev-1] = 2/3*log10(Mo*1e7) - 10.7

    end

    return Mw
end

using StatsBase
using PyPlot

function mfdPlot(Mw)
    
    hist = fit(Histogram, Mw, nbins = 25; closed=:right)

    plot(hist.edges[1][2:end], log10.(hist.weights), ".")
    xlabel("Moment Magnitude")
    ylabel("Log10 of number of earthquakes")
    title("Magnitude frequency distribution")
    show()

end


