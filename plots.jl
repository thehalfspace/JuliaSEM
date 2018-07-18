#################################
#   PLOTS FOR EARTHQUAKE CYCLES
#################################

using PyPlot

# Plot friction parameters
function fricPlot(cca, ccb, FltX)

    plot(cca, FltX/1e3, "k-", label="a")
    plot(ccb, FltX/1e3, "k--", label="b")
    plot(cca-ccb, FltX/1e3, "r-", label="a-b")
    xlabel("Value")
    ylabel("Depth (km)")
    title("Rate and State Friction Parameters")
    legend(loc="upper right")
    ylim([-40, 0])
    show()
end

# Plot shear stress at location as a function of time
function stressPlot(Stress, time_, yr2sec)
    
    loc1 = 3e3  # 3 km depth
    FltID = find(abs.(FltX) .< loc1)[1]

    plot(time_/yr2sec, Stress[FltID, :])
    ylabel("Shear stress at location 1")
    xlabel("time (yr)")
    title("Shear stress as a function of time")
    show()

end


# Plot slip velocity at location as a function of time (same location)
function slipvelPlot(SlipVel, time_, yr2sec)
    
    loc1 = 3e3  # 3 km depth
    FltID = find(abs.(FltX) .< loc1)[1]

    plot(time_/yr2sec, log10.(SlipVel[FltID, :]))
    ylabel("Log of Slip rate at location 1")
    xlabel("time (yr)")
    title("Slip rate as a function of time")
    show()

end


# Plot cumulative slip
function cumSlip(delfsec, delf5yr, FltX)

    indx = find(abs.(FltX) .<= 18e3)[1] 

    plot(delf5yr, FltX/1e3, "b-", linewidth=0.2)
    plot(delfsec, FltX/1e3, "r-", linewidth=0.2)
    xlabel("Slip (m)")
    ylabel("Depth (km)")
    title("Cumulative Slip")
    ylim([-40, 0])

    show()

end
