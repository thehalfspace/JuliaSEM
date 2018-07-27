# Earthquake magnitude Calculation

# Known parameters:
#   SlipVel: Slip velocity at every point, every timestep
#   Slip: Accumulated slip at every point, every timestep

# Compute maximum sliprate and location of maximum sliprate
# at every timestep
function Vfmax(SlipVel, it)
    Vfmax = zeros(it)
    idx = zeros(it)
    for eo = 1:it-1
        Vfmax[eo], idx[eo] = findmax(SlipVel[:,eo])
    end

    return Vfmax, idx
end



# Compute and plot number of seismic events
function events(Vfmax, time_)
    
