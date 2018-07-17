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
    show()
end

# Plot shear stress at location as a function of time
loc1 = 3e3  # 3 km depth
FltID = find(abs.(FltX) .< loc1)[1]


# Plot slip velocity at location as a function of time (same location)

