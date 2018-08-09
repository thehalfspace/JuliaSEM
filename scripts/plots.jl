#################################
#   PLOTS FOR EARTHQUAKE CYCLES
#################################

using PyPlot

# Customize plot
#PyPlot.matplotlib[:rc]("text", usetex = true)
#PyPlot.matplotlib[:rc]("text.latex", preamble = "\\usepackage{amsmath}")
#PyPlot.matplotlib[:rc]("font", size = 18)
#PyPlot.matplotlib[:rc]("axes", labelsize = 15)
#PyPlot.matplotlib[:rc]("axes", titlesize = 15)
#PyPlot.matplotlib[:rc]("xtick", labelsize = 12)
#PyPlot.matplotlib[:rc]("ytick", labelsize = 12)
#PyPlot.matplotlib[:rc]("legend", fontsize = 18)
#PyPlot.matplotlib[:rc]("figure", titlesize = 18)
#PyPlot.matplotlib[:rc]("figure", figsize = (6,4), dpi = 120)

# Plot friction parameters
function fricPlot(cca, ccb, FltX)
    
    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](cca, FltX/1e3, "k--", label="a", lw = 1)
    ax[:plot](ccb, FltX/1e3, "k--", label="b", lw = 1)
    ax[:plot](cca-ccb, FltX/1e3, label="a-b", lw = 1)
    ax[:set_xlabel]("Value")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Rate and State friction parameters")
    ax[:legend](loc="upper right")
    ax[:set_ylim]([-24, 0])
    show()

    figname = string(dir, "/plots", name, "/fric.png")
    fig[:savefig](figname, dpi = 300)

end

# Plot shear stress at location as a function of time
function stressPlot(Stress, time_, FltX, yr2sec, loc1 = 8e3)
    
    FltID = find(abs.(FltX) .<= loc1)[1]
    
    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_/yr2sec, Stress[FltID, :], lw = 1)
    ax[:set_xlabel]("Time (years)")
    ax[:set_ylabel]("Shear stress (MPa)")
    ax[:set_title](string("Shear stress at ", loc1/1e3, "km depth"))   
    show()

    figname = string(dir, "/plots", name, "/shear.png")
    fig[:savefig](figname, dpi = 300)

end


# Plot slip velocity at location as a function of time (same location)
function slipvelPlot(SlipVel, time_, FltX, yr2sec, loc1 = 8e3)
    
    FltID = find(abs.(FltX) .<= loc1)[1]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_/yr2sec, SlipVel[FltID, :], lw = 1)
    ax[:set_xlabel]("Time (years)")
    ax[:set_ylabel]("Slip rate (m/s)")
    ax[:set_title](string("Slip rate at ", loc1/1e3, "km depth"))
    ax[:set_yscale]("log")
    show()

    figname = string(dir, "/plots", name, "/sliprate.png")
    fig[:savefig](figname, dpi = 300)

end


# Plot cumulative slip
function cumSlip(delfsec, delf5yr, FltX)

    indx = find(abs.(FltX) .<= 18e3)[1]

    delfsec2 = delfsec[indx:end, :]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](delf5yr, FltX/1e3, "xkcd:blue", lw=0.5)
    ax[:plot](delfsec2, FltX[indx:end]/1e3, "xkcd:burnt orange", lw=0.5)
    ax[:set_xlabel]("Slip (m)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Cumulative Slip for interseismic and dynamic events")
    ax[:set_ylim]([-24, 0])
    ax[:set_xlim]([0, maximum(delfsec)])
    show()
    
    figname = string(dir, "/plots", name, "/cumslip.png")
    fig[:savefig](figname, dpi = 300)

end
