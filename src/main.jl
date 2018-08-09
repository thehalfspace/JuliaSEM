#########################################
#										
#	SPECTRAL ELEMENT METHOD FOR  		
#	EARTHQUAKE CYCLE SIMULATION			
#										
#	Prithvi Thakur						
#	Adapted from Kaneko et al. (2011)	
#	and J.P. Ampuero's SEMLAB       	
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
include("PCG.jl")
include("dtevol.jl")
include("NRsearch.jl")
include("otherFunctions.jl")


function main(s::space_parameters, tim::time_parameters, 
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
    dt = tim.CFL*dt
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2

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
    #tauAB::Array{Float64} = zeros(s.FltNglob)

    # Initial state variable
    psi::Array{Float64} = tauo./(Seff.*ccb) - eq.fo./ccb - (cca./ccb).*log.(2*v[iFlt]./eq.Vo)
    psi0::Array{Float64} = psi[:]

    # Compute XiLF used in timestep calculation
    XiLf = XiLfFunc(s, tim, eq, muMax, cca, ccb, Seff) 


    # Time related non-constant variables variables
    #slipstart::Int = 0
    #ievb::Int = 0
    #ieva::Int = 0
    ntvsx::Int = 0
    nevne::Int = 0
    isolver::Int = 1
    tvsx::Int64 = 2*tim.yr2sec
    tvsxinc::Int64 = tvsx
    
    # Skip lines 486-490
    # Skip lines 492-507: Outloc1, 2, variables.

    # Display important parameters
    println("Total number of nodes on fault: ", s.FltNglob)
    println("Average node spacing: ", s.LX/(s.FltNglob-1))
    @printf("dt: %1.09f s\n", dt)

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

    while t < tim.Total_time
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
                
                psi1, Vf1 = slrFunc(eq, NFBC, s.FltNglob, psi, psi1, Vf, Vf1, 
                                    IDstate, tau1, tauo, Seff, cca, ccb, dt)

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
            
            dPre .= d[:]
            vPre .= v[:]

            # Update
            d .= d .+ dt.*v .+ (half_dt_sq).*a

            # Prediction
            v .= v .+ half_dt.*a
            a[:] .= 0

            # Internal forces -K*d[t+1] stored in global array 'a'
            # This is different from matlab code; will change if Nel_ETA is not zero
            a = element_computation2(s, iglob, d, H, Ht, W, a)
            a[FltIglobBC] .= 0

            # Absorbing boundaries
            a[iBcL] .= a[iBcL] .- BcLC.*v[iBcL]
            a[iBcT] .= a[iBcT] .- BcTC.*v[iBcT]

            ###### Fault Boundary Condition: Rate and State #############
            FltVfree .= 2*v[iFlt] .+ 2*half_dt*a[iFlt]./M[iFlt]
            Vf .= 2*vPre[iFlt] .+ eq.Vpl


            #for jF = 1:FaultNglob-NFBC
            for j = NFBC: s.FltNglob-1 

                #j = jF - 1 + NFBC
                psi1[j] = IDS(eq.xLf[j], eq.Vo[j], psi[j], dt, Vf[j], 1e-5, IDstate)

                Vf1[j], tau1[j] = NRsearch(eq.fo[j], eq.Vo[j], cca[j], ccb[j],Seff[j],
                                          tauNR[j], tauo[j], psi1[j], FltZ[j], FltVfree[j])
            
                if Vf[j] > 1e10 || isnan(Vf[j]) == 1 || isnan(tau1[j]) == 1
                    #println(FltVfree)
                    println("Fault location = ", j)
                    
                    # Directory to save the simulation results
                    filename = string(dir, "/data", name, "nrfail.jld")

                    @save filename
                    error("NR SEARCH FAILED!")
                    return
                end
                
                psi2[j] = IDS2(eq.xLf[j], eq.Vo[j], psi[j], psi1[j], dt, Vf[j], Vf1[j], IDstate)
                
                # NRsearch 2nd loop
                Vf2[j], tau2[j] = NRsearch(eq.fo[j], eq.Vo[j], cca[j], ccb[j],Seff[j],
                                          tau1[j], tauo[j], psi2[j], FltZ[j], FltVfree[j])

            end
            
            tau .= tau2[:] .- tauo[:]
            tau[iFBC] .= 0
            psi .= psi2[:]
            #KD = a[:]
            a[iFlt] .= a[iFlt] - FltB.*tau
            ########## End of fault boundary condition ############## 


            #RHS = a[:]

            # Solve for a_new
            a[:] .= a./M
            
            # Correction
            v .= v .+ half_dt*a

            v[FltIglobBC] .= 0
            a[FltIglobBC] .= 0

            #### Line 861: Omitting P_Ma
            
            #LHS = M.*a
            #RMS = sqrt.(sum.((RHS - LHS).^2)/length(RHS))./maximum(abs.(RHS))

        end # of isolver if loop
        
        Vfmax = 2*maximum(v[iFlt]) + eq.Vpl


        #----
        # Output variables at different depths for every timestep
        # Omitted the part of code from line 871 - 890, because I 
        # want to output only certain variables each timestep
        #----

        # Output stress, slip, sliprate on fault every certain interval
        if t > tvsx
            ntvsx = ntvsx + 1
            
            delf5yr[:,ntvsx] = 2*d[iFlt] + eq.Vpl*t
            Vf5yr[:,ntvsx] = 2*v[iFlt] + eq.Vpl
            Tau5yr[:,ntvsx] = (tau + tauo)./1e6
            
            tvsx = tvsx +tvsxinc
        end
        
        if Vfmax > eq.Vevne 
            if idelevne == 0
                nevne = nevne + 1
                idelevne = 1
                tevneb = t
                tevne = tim.tevneinc

                delfsec[:,nevne] = 2*d[iFlt] + eq.Vpl*t
                Vfsec[:,nevne] = 2*v[iFlt] + eq.Vpl
                Tausec[:,nevne] = (tau + tauo)./1e6
            end

            if idelevne == 1 && (t - tevneb) > tevne
                nevne = nevne + 1
                
                delfsec[:,nevne] = 2*d[iFlt] + eq.Vpl*t
                Vfsec[:,nevne] = 2*v[iFlt] + eq.Vpl
                Tausec[:,nevne] = (tau + tauo)./1e6

                tevne = tevne + tim.tevneinc
            end

        else
            idelevne = 0
        end

        #-----
        # Output stress and slip before and after events
        # Omitting lines 920-934
        #-----

        # Output timestep info on screen
        if mod(it,500) == 0
            @printf("\nTime (yr) = %1.5g", t/tim.yr2sec)
        end
        
        # Determine quasi-static or dynamic regime based on max-slip velocity
        if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
            isolver = 1
        else
            isolver = 2
        end


        # Some variables for each timestep
        Stress[:,it] = (tau + tauo)./1e6
        SlipVel[:,it] = 2*v[iFlt] + eq.Vpl
        Slip[:,it] = 2*d[iFlt] + eq.Vpl*t
        
        # Compute next timestep dt
        dt = dtevol(tim, dt , dtmin, XiLf, s.FltNglob, NFBC, SlipVel[:,it], isolver)

    end # end of time loop
    
    # Remove zeros from preallocated vectors
    time_ = time_[1:it]

    delfsec = delfsec[:, 1:nevne]
    Vfsec = Vfsec[:,1:nevne]
    Tausec = Tausec[:,1:nevne]

    delf5yr = delf5yr[:,1:ntvsx]
    Vf5yr = Vf5yr[:,1:ntvsx]
    Tau5yr = Tau5yr[:,1:ntvsx]

    Stress = Stress[:,1:it]
    SlipVel = SlipVel[:,1:it]
    Slip = Slip[:,1:it]

    return FltX, delf5yr, delfsec, Stress, SlipVel, Slip, time_, cca, ccb 
end
