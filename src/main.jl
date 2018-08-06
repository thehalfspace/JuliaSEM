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
include("PCG.jl")
include("dtevol.jl")
include("NRsearch.jl")
include("otherFunctions.jl")


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
    Vf0::Array{Float64} = zeros(length(iFlt))
    FltVfree::Array{Float64} = zeros(length(iFlt))
    psi1::Array{Float64} = zeros(s.FltNglob)
    psi2::Array{Float64} = zeros(s.FltNglob)
    tau1::Array{Float64} = zeros(s.FltNglob)
    tau2::Array{Float64} = zeros(s.FltNglob)
    tau3::Array{Float64} = zeros(s.FltNglob)
    tauNR::Array{Float64} = zeros(s.FltNglob)
    tauAB::Array{Float64} = zeros(s.FltNglob)

    # Initial state variable
    psi::Array{Float64} = tauo./(Seff.*ccb) - eq.fo./ccb - (cca./ccb).*log.(2*v[iFlt]./eq.Vo)
    psi0::Array{Float64} = psi[:]

    # Compute XiLF used in timestep calculation
    XiLf = XiLfFunc(s, t, eq, muMax, cca, ccb, Seff) 


    # Time related non-constant variables variables
    slipstart::Int = 0
    ievb::Int = 0
    ieva::Int = 0
    ntvsx::Int = 0
    nevne::Int = 0
    isolver::Int = 1
    
    # Skip lines 486-490
    # Skip lines 492-507: Outloc1, 2, variables.

    # Display important parameters
    println("Total number of nodes on fault: ", s.FltNglob)
    println("Average node spacing: ", s.LX/(s.FltNglob-1))
    @printf("dt: %1.09f s", dt)

    # Find nodes that do not belong to the fault
    FltNI = deleteat!(collect(1:nglob), iFlt)

    # Some more initializations
    r::Array{Float64} = zeros(nglob)
    beta_::Array{Float64} = zeros(nglob)
    alpha_::Array{Float64} = zeros(nglob)
    p::Array{Float64} = zeros(nglob)

    F::Array{Float64} = zeros(nglob)
    dPre::Array{Float64} = zeros(nglob)
    vPre::Array{Float64} = zeros(nglob)
    dd::Array{Float64} = zeros(nglob)
    dacum::Array{Float64} = zeros(nglob)
    
    dnew::Array{Float64} = zeros(length(FltNI))

    
    # Compute diagonal of K
    diagKnew = KdiagFunc(s, iglob, W, H, Ht, FltNI) 

    v[:] = v[:] - 0.5*eq.Vpl
    Vf::Array{Float64} = 2*v[iFlt]
    iFBC::Array{Int64} = find(abs.(FltX) .> 24e3)
    NFBC::Int64 = length(iFBC)
    Vf[iFBC] = 0

    # Fault boundary: indices where fault within 24 km
    fbc = reshape(iglob[:,1,:], length(iglob[:,1,:]),1)
    idx = find(fbc .== find(x .== -24e3)[1] - 1)[1]
    FltIglobBC::Array{Int64} = fbc[1:idx]

    v[FltIglobBC] = 0


    # Preallocate variables with unknown size
    time_ = zeros(1e6)

    delfsec::Array{Float64} = zeros(s.FltNglob, 1e5)
    Vfsec::Array{Float64} = zeros(s.FltNglob, 1e5)
    Tausec::Array{Float64} = zeros(s.FltNglob, 1e5)

    delf5yr::Array{Float64} = zeros(s.FltNglob, 1e4)
    Vf5yr::Array{Float64} = zeros(s.FltNglob, 1e4)
    Tau5yr::Array{Float64} = zeros(s.FltNglob, 1e4)

    Stress::Array{Float64} = zeros(s.FltNglob, 1e6)
    SlipVel::Array{Float64} = zeros(s.FltNglob, 1e6)
    Slip::Array{Float64} = zeros(s.FltNglob, 1e6)


    #....................
    # Start of time loop
    #....................
    it = 0
    t = 0
    IDstate = 2

    #while t < t.Total_time
    while it<10
        it = it + 1
        t = t + dt

        time_[it] = t 


        if isolver == 1

            vPre .= v[:]
            dPre .= d[:]

            Vf0 .= 2*v[iFlt] + eq.Vpl
            Vf  .= Vf0[:]

            for p1 = 1:2
                
                # Compute the forcing term
                F[:] .= 0
                F[iFlt] .= dPre[iFlt] .+ v[iFlt]*dt

                # Assign previous displacement field as initial guess
                dnew .= d[FltNI]

                # Solve d = K^-1F by PCG
                dnew = PCG(s, diagKnew, dnew, F, iFlt, FltNI,
                              H, Ht, iglob, nglob, W)
                
                # update displacement on the medium
                d[FltNI] .= dnew

                # make d = F on the fault
                d[iFlt] .= F[iFlt]

                # Compute on-fault stress
                a[:] .= 0

                a = element_computation(s, iglob, d, H, Ht, W, a)

                a[FltIglobBC] .= 0
                tau1 .= -a[iFlt]./FltB
                
                # Compute slip-rates on-fault
                for j = NFBC: s.FltNglob-1 

                    psi1[j] = IDS(psi[j], dt, eq.Vo[j], eq.xLf[j], Vf[j], 1e-6, IDstate)

                    tauAB[j] = tau1[j] + tauo[j]
                    fa = tauAB[j]/(Seff[j]*cca[j])
                    help = -(eq.fo[j] + ccb[j]*psi1[j])/cca[j]
                    help1 = exp(help + fa)
                    help2 = exp(help - fa)
                    Vf1[j] = eq.Vo[j]*(help1 - help2) 
                end
                
                Vf1[iFBC] .= eq.Vpl
                Vf .= (Vf0 + Vf1)/2
                v[iFlt] .= 0.5*(Vf - eq.Vpl)

            end

            psi .= psi1[:]
            tau .= tau1[:]
            tau[iFBC] .= 0
            Vf1[iFBC] .= eq.Vpl

            v[iFlt] .= 0.5*(Vf1 - eq.Vpl)
            v[FltNI] .= (d[FltNI] - dPre[FltNI])/dt

            #RHS = a[:]
            #RHS[iFlt] = RHS[iFlt] - FltB.*tau
            #RMS = sqrt(sum(RHS.^2)/length(RHS))./maximum(abs.(RHS))
            
            # Line 731: P_MA: Omitted
            a[:] .= 0
            d[FltIglobBC] .= 0
            v[FltIglobBC] .= 0

            
            # If isolver != 1, or max slip rate is < 10^-2 m/s
        else
            a = 1
        end

    end

    return FltIglobBC
end

s = space_parameters()
t = time_parameters()
m = medium_properties()
eq = earthquake_parameters()
FltIglobBC = main(s,t,m,eq);
