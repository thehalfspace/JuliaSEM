#################################
#   PLOTS FOR EARTHQUAKE CYCLES
#################################

using PyPlot

# Customize plot
#PyPlot.matplotlib[:rc]("text", usetex = true)
#PyPlot.matplotlib[:rc]("text.latex", preamble = "\\usepackage{amsmath}")
#  PyPlot.matplotlib[:rc]("font", size = 24)
#  PyPlot.matplotlib[:rc]("axes", labelsize = 21)
#  PyPlot.matplotlib[:rc]("axes", titlesize = 21)
#PyPlot.matplotlib[:rc]("xtick", labelsize = 12)
#PyPlot.matplotlib[:rc]("ytick", labelsize = 12)
#  PyPlot.matplotlib[:rc]("legend", fontsize = 15)
#  PyPlot.matplotlib[:rc]("figure", titlesize = 24)
#PyPlot.matplotlib[:rc]("figure", figsize = (6,4), dpi = 120)


#  function 3DSliprates(S, SlipVel, time_)

    #  # depth and time intervals
    #  sdx = 10; sdt = 100

    #  fig = PyPlot.figure(figsize=(6, 4.5), dpi=120)

#  end

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

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "friction.png")
    fig[:savefig](figname, dpi = 300)

end

# Plot shear stress at location as a function of time
function stressPlot(Stress, time_, FltX, yr2sec, loc1 = 8e3)
    
    FltID = findall(abs.(FltX) .<= loc1)[1]
    
    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_/yr2sec, Stress[FltID, :], lw = 1)
    ax[:set_xlabel]("Time (years)")
    ax[:set_ylabel]("Shear stress (MPa)")
    ax[:set_title](string("Shear stress at ", loc1/1e3, "km depth"))   
    show()

    #  figname = string(dir, "/plots", name, "/shear.png")
    #  fig[:savefig](figname, dpi = 300)

end


# Plot slip velocity at location as a function of time (same location)
function slipvelPlot(SlipVel, time_, FltX, yr2sec, loc1 = 8e3)
    
    FltID = findall(abs.(FltX) .<= loc1)[1]

    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_/yr2sec, SlipVel[FltID, :], lw = 1)
    ax[:set_xlabel]("Time (years)")
    ax[:set_ylabel]("Slip rate (m/s)")
    ax[:set_title](string("Slip rate at ", loc1/1e3, "km depth"))
    ax[:set_yscale]("log")
    show()

    #  figname = string(dir, "/plots", name, "/sliprate.png")
    #  fig[:savefig](figname, dpi = 300)

end

# Plot Vfmax
function VfmaxPlot(SlipVel, time_, yr2sec)

    Vfmax = maximum(SlipVel, dims = 1)[:]
    
    fig = PyPlot.figure(figsize=(6,4.5), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_./yr2sec, Vfmax, lw = 1)
    ax[:set_xlabel]("Time (years)")
    ax[:set_ylabel]("Max. Slip rate (m/s)")
    ax[:set_title]("Max. slip rate on fault")
    ax[:set_yscale]("log")
    show()

    figname = string(path, "Vfmax.png")
    fig[:savefig](figname, dpi = 300)
    #  figname = string(dir, "/plots", name, "/Vfmax.png")
    #  fig[:savefig](figname, dpi = 300)
end


# Plot cumulative slip
function cumSlipPlot(delfsec, delf5yr, FltX)
    
    FZ = -8

    indx = findall(abs.(FltX) .<= 18e3)[1]

    delfsec2 = delfsec[indx:end, :]

    fig = PyPlot.figure()
    ax = fig[:add_subplot](111)
    
    # Shade the fault zone region
    x_shade = LinRange(5,25,25)
    y1 = repeat([0],25)
    y2 = repeat([FZ],25)
    y3 = repeat([-24],25)

    ax[:plot](delf5yr, FltX/1e3, color="royalblue", lw=1, alpha=1.0)
    #  ax[:plot](delfsec2, FltX[indx:end]/1e3, "--", color="chocolate", lw=1, alpha=1.0)
    #  ax[:fill_between](x_shade, y2, y1, color="chocolate", alpha=0.3)
    #  ax[:fill_between](x_shade, y3, y1, color="chocolate", alpha=0.3)
    ax[:set_xlabel]("Accumulated Slip (m)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Cumulative Slip History")
    ax[:set_ylim]([-24, 0])
    ax[:set_xlim]([5,10])  #[0,maximum(delf5yr)])
    show()
    
    figname = string(path, "cumslip_demo2.pdf")
    fig[:savefig](figname, dpi = 300)

end

function stress_slip(O, FltX, start_id, end_id, depth=8e3)
    id = findmin(FltX .< -depth)[2]
    
    fig = PyPlot.figure(figsize=(12,8))
    ax = fig[:add_subplot](111)

    ax[:plot](O.Slip[id, start_id:end_id], O.Stress[id, start_id:end_id], "k.")
    ax[:set_xlabel]("Cumulative Slip (m)")
    ax[:set_ylabel]("Shear Stress (MPa)")
    ax[:set_title]("Slip vs. frictional stress")
    show()

    figname = string(path, "sliptraction.png")
    fig[:savefig](figname, dpi=300)

end

