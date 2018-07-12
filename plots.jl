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
