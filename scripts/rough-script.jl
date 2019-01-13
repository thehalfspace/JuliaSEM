#################################
# SOME ROUGH SCRIPTS
# FOR TRYING OUT STUFF
#################################

using StatsBase
using PyPlot

# Get index for each event
function event_indx(tStart, tEnd, time_)
    start_indx = zeros(size(tStart))
    end_indx = zeros(size(tEnd))

    for i = 1:length(tStart)
        
        start_indx[i] = findall(time_ .== tStart[i])[1]
        end_indx[i] = findall(time_ .== tEnd[i])[1]
    end

    return start_indx, end_indx
end

# Plot sliprates for each event with depth
function test1(S, time_, start_indx, end_indx, SlipVel)
    start_indx = Int.(start_indx)
    end_indx = Int.(end_indx)
    sv = zeros(size(SlipVel))
    inv = 1       # time interval = 0.1 sec for plotting
    to = time_[start_indx]    # start time
    j = 1
    for i=start_indx:end_indx
        if time_[i] >= to
            sv[:,j] = SlipVel[:,i]
            to = to + inv
            j = j+1
        end
    end
    
    sv = sv[:,1:j]

    fig = PyPlot.figure(figsize=(8,6))
    ax = fig[:add_subplot](111)
    
    ax[:plot](sv, S.FltX/1e3, "k--", label="a", lw = 1)
    ax[:set_xlabel]("Shear Stress (MPa)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Shear stress for one event Mw ~ 7")
    #  ax[:set_xlim]([0, 0.02])
    ax[:set_ylim]([-24, 0])
    show()
    
    figname = string(path, "stressM7.pdf")
    fig[:savefig](figname, dpi = 300)
end


# Trying out plotting vertices with imshow
function patch(x, y, v)
    x .= x./1e3
    y .= y./1e3
    LX = maximum(x)
    LY = maximum(y)
    fig = PyPlot.figure(figsize=(8,6))
    ax = fig[:add_subplot](111)
    
    imm = ax[:imshow](v, vmin = minimum(abs.(v)), vmax = maximum(abs.(v)), 
               extent = [0,1,0,1], interpolation="bilinear")
    PyPlot.colorbar(imm)
    show()

end
