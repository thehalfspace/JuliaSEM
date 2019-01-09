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

# Plot shear stress at the start of each event
start_indx, end_indx = event_indx(tStart, tEnd, O.time_)
si = Int.(start_indx)
function stress_event(si, Stress,FltX)
    
    fig = PyPlot.figure(figsize=(12,8), dpi = 120)
    ax = fig[:add_subplot](111)

    ax[:plot](Stress[:,si], FltX./1e3, "k--", lw = 1)
    ax[:set_xlabel]("Shear Stress (MPa)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("Shear stress at the start of each event")
    ax[:set_ylim]([-24, 0])
    show()

    #  figname = string(dir, "/plots", name, "/fric.png")
    figname = string(path, "stress_event.pdf")
    fig[:savefig](figname, dpi = 300)

end


# Animate the rupture events: stresses and sliprates
using PyCall
pygui(true)
@pyimport matplotlib.animation as animation

function stress_movie2(Stress, start_indx, end_indx, evno, FltX)
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
    
    function animate(i)
        line[:set_data](data[:,i], FltX)
        return(line,"")
    end

    ani = animation.FuncAnimation(fig, animate, blit=false, interval=10, frames=200, repeat=false)

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
    tStart = Int(tStart)
    tEnd = Int(tEnd)
    sv = SlipVel[:,tStart:n:tEnd]
    
    fig = PyPlot.figure(figsize=(8,6))
    ax = fig[:add_subplot](111)
    
    ax[:plot](sv, S.FltX/1e3, "k--", label="a", lw = 1)
    ax[:set_xlabel]("Shear Stress (MPa)")
    ax[:set_ylabel]("Depth (km)")
    ax[:set_title]("1.5 km wide fault zone. Mw = 4")
    #  ax[:set_xlim]([0, 0.02])
    ax[:set_ylim]([-24, 0])
    show()
    
    figname = string(path, "shearstressM4_1500.png")
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
