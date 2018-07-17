#################################
#   PLOTS FOR EARTHQUAKE CYCLES
#################################

using PyPlot

# Plot friction parameters
function fricPlot(cca, ccb, FltX)

    plot(cca, FltX/1e3, "k-")
    plot(ccb, FltX/1e3, "k--")
    plot(cca-ccb, FltX/1e3, "r-")
    plot(zeros(cca), FltX/1e3, "k-.")
    show()
end


# Plot shear stress at location as a function of time
loc1 = 3e3  # 3 km depth
FltID = find(abs.(FltX) .< loc1)[1]

shearPlot = plot(x = time_./yr2sec, y = Stress[FltID,:], Geom.line,
                Guide.xlabel("Time (yr)"),
                Guide.ylabel("Shear stress at location 1"),
                Guide.title("Shear Stress as a function of time"),
                Coord.Cartesian(ymin=10, ymax = 60, xmax = 1000))


# Plot slip velocity at location as a function of time (same location)
slipVelPlot = plot(x = time_./yr2sec, y = SlipVel[FltID,:], Geom.line,
                Guide.xlabel("Time (yr)"),
                Guide.ylabel("Shear stress at location 1"),
                Guide.title("Shear Stress as a function of time"),
                Coord.Cartesian(ymin=0, ymax = 5, xmax = 1000))
