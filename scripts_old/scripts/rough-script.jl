#################################
# SOME ROUGH SCRIPTS
# FOR TRYING OUT STUFF
#################################

using StatsBase
using Statistics

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

# Plot shear stress at the start of each event
start_indx, end_indx = event_indx(tStart, tEnd, O.time_)
function stress_avg(Stress,time_, tStart)
   
    str = mean(Stress, dims=1)[:]
    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_, str, "k-", lw = 1)
    for i=1:length(tStart)
        ax[:axvline](tStart[i], linestyle="--", color="black", alpha=0.4)
    end

    ax[:set_xlabel]("Time (yr)")
    ax[:set_ylabel]("Average Shear Stress on Fault (MPa)")
    ax[:set_title]("Average Shear Stress in Time")
    ax[:set_xlim]([471.76,471.765])
    ax[:set_ylim]([24,25])
    show()

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "stress_avg_zoom.pdf")
    fig[:savefig](figname, dpi = 300)

end

si = Int.(start_indx)
function stress_event(a1, a2, Stress,FltX)
   
    str = Stress[:, a1:60:a2]
    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](str, FltX./1e3, "-", lw = 1)
    ax[:set_xlabel]("Slip Velocity (m/s)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Shear stress at the start of each event")
    ax[:set_ylim]([-24, 0])
    show()

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "sliprate_event.pdf")
    #  fig[:savefig](figname, dpi = 300)

end

function stress_max(Stress,FltX, time_)
   
    str, strid = findmax(Stress, dims=1)
    strid = strid[:]
    y_axis = zeros(length(str))
    for i in CartesianIndices(strid)
        y_axis[i] = FltX[strid[i][1]]
    end

    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](time_, y_axis./1e3, "ko", markersize=3, lw = 1)
    #  ax[:plot](collect(1:1000), -8*ones(1000), "k--", label="Fault Zone Depth")
    ax[:set_xlabel]("Time (yr)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Shear stress peaks through time")
    ax[:set_ylim]([-24, 0])
    #  ax[:legend](loc="upper right")
    show()

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "stress_peak.pdf")
    fig[:savefig](figname, dpi = 300)

end

# Animate the rupture events: stresses and sliprates
using PyCall
pygui(true)
@pyimport matplotlib.animation as anim

function stress_movie2(Stress, start_indx, end_indx, evno, FltX)
    indx1 = Int(start_indx[evno])
    indx2 = Int(end_indx[evno])

    data = Stress[:,indx1:indx1+10]

    fig = PyPlot.figure(figsize=(8,6))
    ax = fig[:add_subplot](111)
    ax[:set_xlabel]("Shear Stress (MPa)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Shear Stress Evolution for one earthquake")
    ax[:set_xlim](0,40)
    ax[:set_ylim](-24,0)
    
    line = ax[:plot]([], [], "k--.")[1]
    
    function init()
        line[:set_data]([],[])
        return line
    end

    function animate(i)
        line[:set_data](data[:,i], FltX)
        return line
    end

    ani = anim.FuncAnimation(fig, animate, init_func=init, interval=10, frames=200)

    show() 

end

function stress_movie(Stress, start_indx, end_indx, evno, FltX)
    indx1 = Int(start_indx[evno])
    indx2 = Int(end_indx[evno])

    data = Stress[:,indx1:indx2]

    fig = PyPlot.figure(figsize=(8,6))
    ax = fig[:add_subplot](111)
    ax[:set_xlabel]("Shear Stress (MPa)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Shear Stress Evolution for one earthquake")
    ax[:set_xlim](0,40)
    ax[:set_ylim](-24,0)
    
    line = ax[:plot]([], [], "k--.")
    function simData()
        function it()
            for indx=indx1:indx2
                produce(data[:,indx],FltX)
            end
        end
        Task(it)
    end
    
    function animate()
        task = simData()

        function points(frame_number)
            dat, fltx = consume(task)
            line[:set_data](dat, fltx)
            return(line,"")
        end
        points
    end

    ani = animation.FuncAnimation(fig, animate(), blit=false, interval=10, frames=200, repeat=false)

    show() 

end

# Plot sliprates for each event with depth
function test1(S, tStart, tEnd, SlipVel, n)
    tStart = Int.(tStart)
    tEnd = Int.(tEnd)
    sv = SlipVel[:,tStart:n:tEnd]
    
    fig = PyPlot.figure()
    ax = fig[:add_subplot](111)
    
    ax[:plot](sv, S.FltX/1e3, "-", label="a", lw = 1)
    ax[:set_xlabel]("Co-seismic Sliprate")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Sliprate")
    ax[:set_xlim]([0, 2])
    #  ax[:set_xscale]("log")
    ax[:set_ylim]([-24, 0])
    show()
    
    figname = string(path, "shearstressM4_1500.png")
    fig[:savefig](figname, dpi = 300)
end
