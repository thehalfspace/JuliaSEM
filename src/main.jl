#########################################
#										
#	SPECTRAL ELEMENT METHOD FOR  		
#	EARTHQUAKE CYCLE SIMULATION			
#										
#	Prithvi Thakur						
#	Adapted from Kaneko et al. (2011)	
#	and J.P. Ampuero's SEMLAB       	
#########################################

using Printf
using LinearAlgebra
using JLD2

#.................................
# Include external function files
#.................................
include("parameters/defaultParameters.jl")	    #	Set Parameters
include("GetGLL.jl")		    #	Polynomial interpolation
include("Meshbox.jl")		    # 	Build 2D mesh
include("Assemble.jl")          #   Assemble mass and stiffness matrix
include("BoundaryMatrix.jl")    #	Boundary matrices
include("FindNearestNode.jl")   #	Nearest node for output
include("setup.jl")             #   Setup the constants for simulation
include("initialConditions/defaultInitialConditions.jl")
include("PCG.jl")               # Preconditioned conjugate gradient to invert matrix
include("dtevol.jl")            # compute the next timestep
include("NRsearch.jl")          # Newton-rhapson search method to find roots
include("otherFunctions.jl")    # some other functions to solve for friction

struct results

    Stress::Array{Float64,2}
    SlipVel::Array{Float64,2}
    Slip::Array{Float64,2}
    time_::Array{Float64}
end

function main(P::parameters, S::input_variables)

    Ht = S.H'
    wgll2 = S.wgll*S.wgll';
    
    # Time solver variables
    dt = S.dt0
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2

    # dt modified slightly for damping
    if P.ETA != 0
	    dt = dt/sqrt(1 + 2*P.ETA)
    end

    # Initialize kinematic field: global arrays
    global d = zeros(S.nglob)
    global v = zeros(S.nglob)
    v .= 0.5e-3
    global a = zeros(S.nglob)

    #.....................................
    # Stresses and time related variables
    #.....................................
    tau::Array{Float64} = zeros(P.FltNglob)
    FaultC::Array{Float64} = zeros(P.FltNglob)
    Vf1::Array{Float64}  = zeros(P.FltNglob)
    Vf2::Array{Float64} = zeros(P.FltNglob)
    Vf0::Array{Float64} = zeros(length(S.iFlt))
    FltVfree::Array{Float64} = zeros(length(S.iFlt))
    psi::Array{Float64} = zeros(P.FltNglob)
    psi0::Array{Float64} = zeros(P.FltNglob)
    psi1::Array{Float64} = zeros(P.FltNglob)
    psi2::Array{Float64} = zeros(P.FltNglob)
    tau1::Array{Float64} = zeros(P.FltNglob)
    tau2::Array{Float64} = zeros(P.FltNglob)
    tau3::Array{Float64} = zeros(P.FltNglob)
    tauNR::Array{Float64} = zeros(P.FltNglob)
    #tauAB::Array{Float64} = zeros(P.FltNglob)

    # Initial state variable
    psi = S.tauo./(S.Seff.*S.ccb) - P.fo./S.ccb - (S.cca./S.ccb).*log.(2*v[S.iFlt]./P.Vo)
    psi0 .= psi[:]

    # Time related non-constant variables variables
    #slipstart::Int = 0
    #ievb::Int = 0
    #ieva::Int = 0
    ntvsx::Int = 0
    nevne::Int = 0
    isolver::Int = 1
    tvsx::Int64 = 2*P.yr2sec
    tvsxinc::Int64 = tvsx
    
    # Skip lines 486-490
    # Skip lines 492-507: Outloc1, 2, variables.

    # Display important parameters
    println("Total number of nodes on fault: ", P.FltNglob)
    println("Average node spacing: ", P.LX/(P.FltNglob-1))
    @printf("dt: %1.09f s\n", dt)

    # Some more initializations
    r::Array{Float64} = zeros(S.nglob)
    beta_::Array{Float64} = zeros(S.nglob)
    alpha_::Array{Float64} = zeros(S.nglob)
    #p::Array{Float64} = zeros(nglob)

    F::Array{Float64} = zeros(S.nglob)
    dPre::Array{Float64} = zeros(S.nglob)
    vPre::Array{Float64} = zeros(S.nglob)
    dd::Array{Float64} = zeros(S.nglob)
    dacum::Array{Float64} = zeros(S.nglob)
    dnew::Array{Float64} = zeros(length(S.FltNI))

    # Iterators
    idelevne = 3
    tevneb = 0
    tevne = 0

    v = v[:] .- 0.5*P.Vpl
    Vf::Array{Float64} = 2*v[S.iFlt]
    iFBC::Array{Int64} = findall(abs.(S.FltX) .> 24e3)
    NFBC::Int64 = length(iFBC)
    Vf[iFBC] .= 0


    v[S.FltIglobBC] .= 0

    # Preallocate variables with unknown size
    time_ = zeros(1000000)

    delfsec::Array{Float64} = zeros(P.FltNglob, 100000)
    #Vfsec::Array{Float64} = zeros(P.FltNglob, 100000)
    #Tausec::Array{Float64} = zeros(P.FltNglob, 100000)

    delf5yr::Array{Float64} = zeros(P.FltNglob, 10000)
    #Vf5yr::Array{Float64} = zeros(P.FltNglob, 10000)
    #Tau5yr::Array{Float64} = zeros(P.FltNglob, 10000)

    Stress::Array{Float64} = zeros(P.FltNglob, 1000000)
    SlipVel::Array{Float64} = zeros(P.FltNglob, 1000000)
    Slip::Array{Float64} = zeros(P.FltNglob, 1000000)


    #....................
    # Start of time loop
    #....................
    it = 0
    t = 0

    while t < P.Total_time
        it = it + 1
        t = t + dt

        time_[it] = t 


        if isolver == 1

            vPre .= v
            dPre .= d

            Vf0 .= 2*v[S.iFlt] .+ P.Vpl
            Vf  .= Vf0

            @inbounds for p1 = 1:2
                
                # Compute the forcing term
                F .= 0
                F[S.iFlt] .= dPre[S.iFlt] .+ v[S.iFlt]*dt

                # Assign previous displacement field as initial guess
                dnew .= d[S.FltNI]

                # Solve d = K^-1F by PCG
                #println("\nPCG:")
                dnew = PCG(P, S.diagKnew, dnew, F, S.iFlt, S.FltNI,
                              S.H, Ht, S.iglob, S.nglob, S.W)
                
                # update displacement on the medium
                d[S.FltNI] .= dnew

                # make d = F on the fault
                d[S.iFlt] .= F[S.iFlt]

                # Compute on-fault stress
                a .= 0

                # Compute forcing (acceleration) for each element
                #println("\nElement Computation:")
                a = element_computation(P, S.iglob, d, S.H, Ht, S.W, a)

                a[S.FltIglobBC] .= 0
                tau1 .= -a[S.iFlt]./S.FltB
                
                # Function to calculate sliprate
               # println("\nSlr function:")
                psi1, Vf1 = slrFunc(P, NFBC, P.FltNglob, psi, psi1, Vf, Vf1, 
                                    P.IDstate, tau1, S.tauo, S.Seff, S.cca, S.ccb, dt)

                Vf1[iFBC] .= P.Vpl
                Vf .= (Vf0 + Vf1)/2
                v[S.iFlt] .= 0.5*(Vf .- P.Vpl)

            end

            psi .= psi1[:]
            tau .= tau1[:]
            tau[iFBC] .= 0
            Vf1[iFBC] .= P.Vpl

            v[S.iFlt] .= 0.5*(Vf1 .- P.Vpl)
            v[S.FltNI] .= (d[S.FltNI] .- dPre[S.FltNI])/dt

            #RHS = a[:]
            #RHS[iFlt] = RHS[iFlt] - FltB.*tau
            #RMS = sqrt(sum(RHS.^2)/length(RHS))./maximum(abs.(RHS))
            
            # Line 731: P_MA: Omitted
            a .= 0
            d[S.FltIglobBC] .= 0
            v[S.FltIglobBC] .= 0

            
            # If isolver != 1, or max slip rate is < 10^-2 m/s
        else
            
            dPre .= d
            vPre .= v

            # Update
            d .= d .+ dt.*v .+ (half_dt_sq).*a

            # Prediction
            v .= v .+ half_dt.*a
            a .= 0

            # Internal forces -K*d[t+1] stored in global array 'a'
            # This is different from matlab code; will change if Nel_ETA is not zero
            a = element_computation2(P, S.iglob, d, S.H, Ht, S.W, a)
            a[S.FltIglobBC] .= 0

            # Absorbing boundaries
            a[S.iBcL] .= a[S.iBcL] .- S.BcLC.*v[S.iBcL]
            a[S.iBcT] .= a[S.iBcT] .- S.BcTC.*v[S.iBcT]

            ###### Fault Boundary Condition: Rate and State #############
            FltVfree .= 2*v[S.iFlt] .+ 2*half_dt*a[S.iFlt]./S.M[S.iFlt]
            Vf .= 2*vPre[S.iFlt] .+ P.Vpl


            #for jF = 1:FaultNglob-NFBC
            @inbounds for j = NFBC: P.FltNglob-1 

                #j = jF - 1 + NFBC
                psi1[j] = IDS(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-5, P.IDstate)

                Vf1[j], tau1[j] = NRsearch(P.fo[j], P.Vo[j], S.cca[j], S.ccb[j], S.Seff[j],
                                          tauNR[j], S.tauo[j], psi1[j], S.FltZ[j], FltVfree[j])
            
                if Vf[j] > 1e10 || isnan(Vf[j]) == 1 || isnan(tau1[j]) == 1
                    
                    println("Fault Location = ", j)

                    # Save simulation results
                    filename = string(dir, "/data", name, "nrfail.jld")
                    @save filename 
                    @error("NR SEARCH FAILED!")
                    return
                end
                
                psi2[j] = IDS2(P.xLf[j], P.Vo[j], psi[j], psi1[j], dt, Vf[j], Vf1[j], P.IDstate)
                
                # NRsearch 2nd loop
                Vf2[j], tau2[j] = NRsearch(P.fo[j], P.Vo[j], S.cca[j], S.ccb[j], S.Seff[j],
                                          tau1[j], S.tauo[j], psi2[j], S.FltZ[j], FltVfree[j])

            end
            
            tau .= tau2 .- S.tauo
            tau[iFBC] .= 0
            psi .= psi2
            #KD = a[:]
            a[S.iFlt] .= a[S.iFlt] .- S.FltB.*tau
            ########## End of fault boundary condition ############## 


            #RHS = a[:]

            # Solve for a_new
            a .= a./S.M
            
            # Correction
            v .= v .+ half_dt*a

            v[S.FltIglobBC] .= 0
            a[S.FltIglobBC] .= 0

            #### Line 861: Omitting P_Ma
            
            #LHS = M.*a
            #RMS = sqrt.(sum.((RHS - LHS).^2)/length(RHS))./maximum(abs.(RHS))

        end # of isolver if loop
        
        Vfmax = 2*maximum(v[S.iFlt]) .+ P.Vpl


        #----
        # Output variables at different depths for every timestep
        # Omitted the part of code from line 871 - 890, because I 
        # want to output only certain variables each timestep
        #----

        # Output stress, slip, sliprate on fault every certain interval
        if t > tvsx
            ntvsx = ntvsx + 1
            
            delf5yr[:,ntvsx] = 2*d[S.iFlt] .+ P.Vpl*t
            #Vf5yr[:,ntvsx] = 2*v[S.iFlt] .+ P.Vpl
            #Tau5yr[:,ntvsx] = (tau + tauo)./1e6
            
            tvsx = tvsx +tvsxinc
        end
        
        if Vfmax > P.Vevne 
            if idelevne == 0
                nevne = nevne + 1
                idelevne = 1
                tevneb = t
                tevne = P.tevneinc

                delfsec[:,nevne] = 2*d[S.iFlt] .+ P.Vpl*t
                #Vfsec[:,nevne] = 2*v[S.iFlt] .+ P.Vpl
                #Tausec[:,nevne] = (tau + tauo)./1e6
            end

            if idelevne == 1 && (t - tevneb) > tevne
                nevne = nevne + 1
                
                delfsec[:,nevne] = 2*d[S.iFlt] .+ P.Vpl*t
                #Vfsec[:,nevne] = 2*v[S.iFlt] .+ P.Vpl
                #Tausec[:,nevne] = (tau + tauo)./1e6

                tevne = tevne + P.tevneinc
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
            @printf("\nTime (yr) = %1.5g", t/P.yr2sec)
        end
        
        # Determine quasi-static or dynamic regime based on max-slip velocity
        if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
            isolver = 1
        else
            isolver = 2
        end


        # Some variables for each timestep
        Stress[:,it] = (tau + S.tauo)./1e6
        SlipVel[:,it] = 2*v[S.iFlt] .+ P.Vpl
        Slip[:,it] = 2*d[S.iFlt] .+ P.Vpl*t
        
        # Compute next timestep dt
        dt = dtevol(P, dt , dtmin, S.XiLf, P.FltNglob, NFBC, SlipVel[:,it], isolver)

    end # end of time loop
    
    # Remove zeros from preallocated vectors
    time_ = time_[1:it]
    Stress = Stress[:,1:it]
    SlipVel = SlipVel[:,1:it]
    Slip = Slip[:,1:it]

    return results(Stress, SlipVel, Slip, time_)

end
