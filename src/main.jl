#########################################
#										#
#	SPECTRAL ELEMENT METHOD FOR  		#
#	EARTHQUAKE CYCLE SIMULATION			#
#										#
#	Prithvi Thakur						#
#	Adapted from Kaneko et al. (2011)	#
#	and J.P. Ampuero's SEMLAB       	#
#########################################

#.................................
# Include external function files
#.................................
include("parameters/defaultParameters.jl")	    #	Set Parameters
include("GetGLL.jl")		#	Polynomial interpolation
include("Meshbox.jl")		# 	Build 2D mesh
include("Assemble.jl")
include("BoundaryMatrix.jl")    #	Boundary matrices
include("FindNearestNode.jl")   #	Nearest node
include("initialConditions/defaultInitialConditions.jl")
include("IDState.jl") # state variable computation
include("XiLf.jl")


function main(s::space_parameters, t::time_parameters, 
              m::medium_properties, eq::earthquake_parameters)

    #....................
    # 2D Mesh generation
    #....................
    iglob, x, y = MeshBox(s)
    x = x - s.LX
    nglob = length(x)

    # The derivatives of the Lagrange Polynomials were pre-tabulated 
    # xgll = location of the GLL nodes inside the reference segment [-1,1]
    xgll, wgll, H = GetGLL(s.NGLL)
    Ht = H'
    wgll2 = wgll*wgll';

    #.................
    # Initialization
    #.................

    # For internal forces
    W::Array{Float64,3} = zeros(s.NGLL, s.NGLL, s.Nel)

    # Global Mass Matrix
    M::Array{Float64} = zeros(nglob)

    # Mass+Damping matrix
    MC::Array{Float64} = zeros(nglob)

    # Assemble mass and stiffness matrix
    M, W, dt, muMax = assemble(s,m,iglob,M,W)
 
    # Time solver variables
    dt = t.CFL*dt
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2

    dtmax = 100 * 24 * 60*60		# 100 days

    # dt modified slightly for damping
    if m.ETA != 0
	    dt = dt/sqrt(1 + 2*m.ETA)
    end

    # Initialize kinematic field: global arrays
    global d = zeros(nglob)
    global v = zeros(nglob)
    v[:] = 0.5e-3
    global a = zeros(nglob)

    #......................
    # Boundary conditions : absorbing boundaries on 3 sides, fault boundary on one side
    #......................

    # Left boundary
    BcLC, iBcL = BoundaryMatrix(s, m, wgll, iglob, 'L')

    # Right Boundary = free surface: nothing to do

    # Top Boundary
    BcTC, iBcT = BoundaryMatrix(s, m, wgll, iglob, 'T')

    # Mass matrix at boundaries
    Mq = M[:]
    M[iBcL] .= M[iBcL] .+ half_dt*BcLC
    M[iBcT] .= M[iBcT] .+ half_dt*BcTC


    # Dynamic fault at bottom boundary
    FltB, iFlt = BoundaryMatrix(s, m, wgll, iglob, 'B') 

    FltZ::Array{Float64} = M[iFlt]./FltB /half_dt * 0.5
    FltX::Array{Float64} = x[iFlt]


    #......................
    # Initial Conditions
    #......................
    cca, ccb = fricDepth(s, FltX)   # rate-state friction parameters

    Seff = SeffDepth(s, FltX)       # effective normal stress

    tauo = tauDepth(s, FltX)        # initial shear stress

    # Kelvin-Voigt Viscosity
    Nel_ETA::Int = 0
    if m.ETA !=0
        Nel_ETA = s.NelX
        x1 = 0.5*(1 + xgll')
        eta_taper = exp.(-pi*x1.^2)
        eta = m.ETA*dt*repmat([eta_taper], s.NGLL)

    else
        Nel_ETA = 0
    end


    #.....................................
    # Stresses and time related variables
    #.....................................
    tau::Array{Float64} = zeros(s.FltNglob)
    FaultC::Array{Float64} = zeros(s.FltNglob)
    Vf1::Array{Float64}  = zeros(s.FltNglob)
    Vf2::Array{Float64} = zeros(s.FltNglob)
    Vf::Array{Float64} 	= zeros(s.FltNglob)
    psi1::Array{Float64} = zeros(s.FltNglob)
    psi2::Array{Float64} = zeros(s.FltNglob)
    tau1::Array{Float64} = zeros(s.FltNglob)
    tau2::Array{Float64} = zeros(s.FltNglob)
    tau3::Array{Float64} = zeros(s.FltNglob)
    tauNR::Array{Float64} = zeros(s.FltNglob)
    tauAB::Array{Float64} = zeros(s.FltNglob)

    # Initial state variable
    psi = tauo./(Seff.*ccb) - eq.fo./ccb - (cca./ccb).*log.(2*v[iFlt]./eq.Vo)
    psi0 = psi[:]

    # Compute XiLF used in timestep calculation
    XiLf = XiLfFunc(s, t, muMax) 

    return cca, tauo
    
end

s = space_parameters()
t = time_parameters()
m = medium_properties()
eq = earthquake_parameters()
cca, tauo = main(s,t,m,eq)
