############################################
# COMPUTE STATIC STRESS DROP FOR EACH EVENT
############################################
using PyPlot
#.................................................
# Compute the final Coseismic slip for each event
#.................................................
function static_stress_drop(Stress, SlipVel, time_=zeros(1e6))

    Vfmax = maximum(SlipVel, 1)

    stressdrop::Array{Float64,2} = zeros(SlipVel)
    t_catalog::Array{Float64} = zeros(SlipVel[1,:])
    
    Vthres = 0.01 # event threshold
    slipstart = 0
    it = 1; it2 = 1
    stress_ref = zeros(SlipVel[:,1]) # Reference stress: shear stress before the earthquake

    eq_depth = zeros(Int, length(time_))

    for i = 1:length(SlipVel[1,:])

        # Start of each event
        if Vfmax[i] > 1.01*Vthres && slipstart == 0
            stress_ref = Stress[:,i]
            slipstart = 1
            t_catalog[it2] = time_[i]

            eq_depth[it2] = find(SlipVel[:,i] .== Vfmax[i])[1]

            it2 = it2+1
        end

        # End of each event
        if Vfmax[i] < 0.99*Vthres && slipstart == 1
            stressdrop[:,it] = Stress[:,i] - stress_ref
            slipstart = 0
            it = it + 1
        end
    end

    return stressdrop[:,1:it-1], t_catalog[1:it2-1], eq_depth[1:it2-1]
end


#......................................
# Plot max. stress drop with depth
#......................................
function stressdropPlot(stressdrop, FltX, eq_depth)
 
    x_axis = maximum(stressdrop, 1)

    fig = PyPlot.figure(figsize=(6,4.5), dpi=120)
    ax = fig[:add_subplot](111)
    
    ax[:plot](x_axis, FltX[eq_depth])
    ax[:set_xlabel]("Stress Drop (MPa)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Static Stress Drop for Earthquakes")
    show()

    #figname = string(dir, "/plots", name, "/catalogue.png")
    #fig[:savefig](figname, dpi = 300)
end
