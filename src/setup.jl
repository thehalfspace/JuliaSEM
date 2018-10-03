include("$(@__DIR__)/GetGLL.jl")		    #	Polynomial interpolation
include("$(@__DIR__)/MeshBox.jl")		    # 	Build 2D mesh
include("$(@__DIR__)/Assemble.jl")          #   Assemble mass and stiffness matrix
include("$(@__DIR__)/BoundaryMatrix.jl")    #	Boundary matrices
include("$(@__DIR__)/FindNearestNode.jl")   #	Nearest node for output
include("$(@__DIR__)/initialConditions/defaultInitialConditions.jl")
include("$(@__DIR__)/otherFunctions.jl")    # some other functions to solve for friction

struct input_variables
    
    iglob::Array{Int,3}
    nglob::Int
    x::Array{Float64}
    y::Array{Float64}

    xgll::Array{Float64}
    wgll::Array{Float64}
    H::Array{Float64,2}
    Ht::Array{Float64,2}

    W::Array{Float64,3}
    M::Array{Float64}
    MC::Array{Float64}

    BcLC::Array{Float64}
    iBcL::Array{Int}
    BcTC::Array{Float64}
    iBcT::Array{Int}

    FltB::Array{Float64}
    iFlt::Array{Int64}
    FltZ::Array{Float64}
    FltX::Array{Float64}

    cca::Array{Float64}
    ccb::Array{Float64}
    Seff::Array{Float64}
    tauo::Array{Float64}

    Nel_ETA::Float64
    XiLf::Array{Float64}
    diagKnew::Array{Float64}

    FltIglobBC::Array{Int}
    FltNI::Array{Int}

    dt0::Float64

end

function setup(P::parameters)


    #....................
    # 2D Mesh generation
    #....................
    iglob, x, y = MeshBox(P)
    x = x .- P.LX
    nglob = length(x)

    # The derivatives of the Lagrange Polynomials were pre-tabulated 
    # xgll = location of the GLL nodes inside the reference segment [-1,1]
    xgll, wgll, H = GetGLL(P.NGLL)
    Ht = H'
    wgll2 = wgll*wgll';

    #.................
    # Initialization
    #.................

    # For internal forces
    W::Array{Float64,3} = zeros(P.NGLL, P.NGLL, P.Nel)

    # Global Mass Matrix
    M::Array{Float64} = zeros(nglob)

    # Mass+Damping matrix
    MC::Array{Float64} = zeros(nglob)

    # Assemble mass and stiffness matrix
    M, W, dt, muMax = assemble(P,iglob,M,W)
    
    # Time solver variables
    dt = P.CFL*dt
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2
    
    #......................
    # Boundary conditions : absorbing boundaries on 3 sides, fault boundary on one side
    #......................

    # Left boundary
    BcLC, iBcL = BoundaryMatrix(P, wgll, iglob, 'L')

    # Right Boundary = free surface: nothing to do
    BcRC, iBcR = BoundaryMatrix(P, wgll, iglob, 'R')

    # Top Boundary
    BcTC, iBcT = BoundaryMatrix(P, wgll, iglob, 'T')

    # Mass matrix at boundaries
    Mq = M[:]
    M[iBcL] .= M[iBcL] .+ half_dt*BcLC
    M[iBcT] .= M[iBcT] .+ half_dt*BcTC
    M[iBcR] .= M[iBcR] .+ half_dt*BcRC


    # Dynamic fault at bottom boundary
    FltB, iFlt = BoundaryMatrix(P, wgll, iglob, 'B') 

    FltZ::Array{Float64} = M[iFlt]./FltB /half_dt * 0.5
    FltX::Array{Float64} = x[iFlt]


    #......................
    # Initial Conditions
    #......................
    cca, ccb = fricDepth(FltX)   # rate-state friction parameters
    Seff = SeffDepth(FltX)       # effective normal stress
    tauo = tauDepth(FltX)        # initial shear stress

    # Kelvin-Voigt Viscosity
    Nel_ETA::Int = 0
    if P.ETA !=0
        Nel_ETA = P.NelX
        x1 = 0.5*(1 .+ xgll')
        eta_taper = exp.(-pi*x1.^2)
        eta = P.ETA*dt*repeat([eta_taper], P.NGLL)

    else
        Nel_ETA = 0
    end

    # Compute XiLF used in timestep calculation
    XiLf = XiLfFunc(P, muMax, cca, ccb, Seff) 
 
    # Find nodes that do not belong to the fault
    FltNI = deleteat!(collect(1:nglob), iFlt)
    
    # Compute diagonal of K
    diagKnew = KdiagFunc(P, iglob, W, H, Ht, FltNI) 
    
    # Fault boundary: indices where fault within 24 km
    fbc = reshape(iglob[:,1,:], length(iglob[:,1,:])) 
    idx = findall(fbc .== findall(x .== -24e3)[1] - 1)[1]    
    FltIglobBC::Array{Int} = fbc[1:idx]

    S = input_variables(iglob,nglob, x, y, xgll, wgll, H, Ht, W, M, MC, BcLC ,iBcL ,BcTC ,
                        iBcT, FltB, iFlt ,FltZ ,FltX ,cca ,ccb ,Seff ,tauo ,Nel_ETA ,
                        XiLf, diagKnew, FltIglobBC ,FltNI, dt)

    return S

end
