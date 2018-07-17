#################################
#   PLOTS FOR EARTHQUAKE CYCLES
#################################

using Gadfly

# Plot friction parameters
fricPlot = plot(layer(x = (cca-ccb).*1e3 + 10, y = FltX/1e3,
                      Theme(default_color=colorant"red"), Geom.point),
                layer(x = Seff/1e6, y = FltX/1e3, 
                      Theme(default_color=colorant"deepskyblue"), Geom.point),
                layer(x = tauo/1e6, y = FltX/1e3, 
                      Theme(default_color=colorant"green"), Geom.point),
                Guide.xlabel("Scaled (a-b)/ Stress value"),
                Guide.ylabel("Depth (m)"),
                Guide.title("Rate and state friction/Stress"),
                Coord.Cartesian(ymin=-24))

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
