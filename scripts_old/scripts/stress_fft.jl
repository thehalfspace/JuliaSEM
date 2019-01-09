#################################
# Looking at shear stresses in the
# fourier domain
#################################

using PyPlot
using FFTW
using LinearAlgebra

function stress_fft(Stress, FltX, start_indx, end_indx, evno)
    
    i1 = Int(start_indx[evno])
    i2 = Int(end_indx[evno])
    #   str = Stress[:,i1:i2]
    str = Stress[:,Int.(start_indx)]

    stress_amp = abs.(fft(str,2))

    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](stress_amp, FltX./1e3, "ko--", lw = 1)
    ax[:set_xlabel]("Shear Stress amplitude in fourier domain")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Shear stress for all events")
    ax[:set_ylim]([-24, 0])
    show()

end
