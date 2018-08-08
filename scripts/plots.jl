#################################
#   PLOTS FOR EARTHQUAKE CYCLES
#################################

using PyPlot

# Customize plot
PyPlot.matplotlib[:rc]("text", usetex = true)
PyPlot.matplotlib[:rc]("text.latex", preamble = "\\usepackage{amsmath}")
#PyPlot.matplotlib[:rc]("font", size = 18)
#PyPlot.matplotlib[:rc]("axes", labelsize = 15)
#PyPlot.matplotlib[:rc]("axes", titlesize = 15)
#PyPlot.matplotlib[:rc]("xtick", labelsize = 12)
#PyPlot.matplotlib[:rc]("ytick", labelsize = 12)
#PyPlot.matplotlib[:rc]("legend", fontsize = 18)
#PyPlot.matplotlib[:rc]("figure", titlesize = 18)
PyPlot.matplotlib[:rc]("figure", figsize = (6,4), dpi = 120)

# Plot friction parameters
function fricPlot(cca, ccb, FltX)

    plot(cca, FltX/1e3, "k-", label="a", lw = 1)
    plot(ccb, FltX/1e3, "k--", label="b", lw = 1)
    plot(cca-ccb, FltX/1e3, "r-", label="a-b", lw = 1)
    xlabel("Value")
    ylabel("Depth (km)")
    title("Rate and State Friction Parameters")
    legend(loc="upper right")
    ylim([-24, 0])
    show()
end

# Plot shear stress at location as a function of time
function stressPlot(Stress, time_, FltX, yr2sec)
    
    loc1 = 8e3  # 3 km depth
    FltID = find(abs.(FltX) .<= loc1)[1]

    plot(time_/yr2sec, Stress[FltID, :])
    ylabel("Shear stress at location 1")
    xlabel("time (yr)")
    title("Shear stress as a function of time")
    show()

end


# Plot slip velocity at location as a function of time (same location)
function slipvelPlot(SlipVel, time_, FltX, yr2sec)
    
    loc1 = 8e3  # 3 km depth
    FltID = find(abs.(FltX) .<= loc1)[1]

    plot(time_/yr2sec, log10.(SlipVel[FltID, :]))
    ylabel("Log of Slip rate at location 1")
    xlabel("time (yr)")
    title("Slip rate as a function of time")
    show()

end


# Plot cumulative slip
function cumSlip(delfsec, delf5yr, FltX)

    indx = find(abs.(FltX) .<= 18e3)[1]

    delfsec2 = delfsec[indx:end, :]

    plot(delf5yr, FltX/1e3, "b-", linewidth=1)
    plot(delfsec2, FltX[indx:end]/1e3, "r-", linewidth=1)
    xlabel("Slip (m)")
    ylabel("Depth (km)")
    title("Cumulative Slip")
    ylim([-24, 0])
    xlim([0, maximum(delfsec)])

    show()

end
