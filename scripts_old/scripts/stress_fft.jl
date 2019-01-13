#################################
# Looking at shear stresses in the
# fourier domain
#################################

#  using Gadfly
using FFTW
using LinearAlgebra
function stress_fft(Stress, FltX, start_indx, end_indx, evno, time_)
    

    x1, x2 = findmax(FltX .>= -8e3)

    i1 = Int(start_indx[evno])
    i2 = Int(end_indx[evno])
     str = Stress[x2,:]
    #  str = Stress[:,Int.(start_indx)]

    stress_amp = abs.(fft(str,1))
    
    #  p1 = plot(x=stress_amp, y=FltX./1e3, Geom.line,
         #  Guide.title("Shear Stress for the Start of Each Event"),
         #  Guide.xlabel("Shear Stress Amplitude in Fourier Domain"),
         #  Guide.ylabel("Depth (km)"))

    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_, stress_amp, "ko", lw = 1)
    ax[:set_xlabel]("Time (yr)")
    ax[:set_ylabel]("Shear Stress Fourier Amplitude")
    #  ax[:set_xscale]("log")
    ax[:set_yscale]("log")
    ax[:set_title]("Shear stress fourier amplitude at 8km depth")
    ax[:set_ylim]([-24, 0])
    show()
end
