#################################
# Looking at shear stresses in the
# fourier domain
#################################

#  using Gadfly
using FFTW
using LinearAlgebra
function stress_fft(Stress, FltX, start_indx, end_indx, evno, time_)
    #  x1, x2 = findmax(FltX .>= -8e3)
    i1 = Int(start_indx[evno])
    i2 = Int(end_indx[evno])
     #  str = Stress[:,]
    str = Stress[:,Int.(start_indx)]

    stress_amp = abs.(fft(str[1:end-1,:],1))
    
    #  p1 = plot(x=stress_amp, y=FltX./1e3, Geom.line,
         #  Guide.title("Shear Stress for the Start of Each Event"),
         #  Guide.xlabel("Shear Stress Amplitude in Fourier Domain"),
         #  Guide.ylabel("Depth (km)"))

    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](stress_amp, -1e3 ./FltX[1:end-1], "ko", lw = 1)
    ax[:set_xlabel]("Shear Stress in fourier domain")
    ax[:set_ylabel]("Depth frequency (km^{-1})")
    ax[:set_xscale]("log")
    ax[:set_yscale]("log")
    ax[:set_title]("Shear stress fourier amplitude at the start of each event")
    #  ax[:set_ylim]([-24, 0])
    ax[:invert_yaxis]()
    show()
    figname = string(path, "fourier_shear_stress.pdf")
    fig[:savefig](figname, dpi = 300)
end

function stress_fft2(stressdrops, FltX)
    stress_amp = abs.(fft(stressdrops[361:end-1,:],1))

    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](stress_amp, FltX[361:end-1]./1e3, "ko--", lw = 1)
    ax[:set_xlabel]("Stress drops in fourier domain")
    ax[:set_ylabel]("Depth frequency (km^-1)")
    ax[:set_xscale]("log")
    #  ax[:set_yscale]("log")
    ax[:set_title]("Stress drop fourier amplitude for each event")
    #  ax[:set_ylim]([-24, 0])
    #  ax[:invert_yaxis]()
    show()
    figname = string(path, "fourier_stress_drop.pdf")
    fig[:savefig](figname, dpi = 300)
end
