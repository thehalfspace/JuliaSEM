###############################################################################
#										
#	SPECTRAL ELEMENT METHOD FOR EARTHQUAKE CYCLE SIMULATION			
#	
#   Written in: Julia 1.0
#
#	Created: 06/20/2018
#   Author: Prithvi Thakur (Original code by Kaneko et al.)
#
#	Adapted from Kaneko et al. (2011)	
#	and J.P. Ampuero's SEMLAB       	
#
#   CHANGELOG:
#       * 08-26-2018: Using distributed for loop in PCG and NRsearch
#       * 08-24-2018: Create a separate function for NRsearch loop: FBC()
#       * 08-20-2018: Use JLD2 to store data instead of JLD
#       * 08-14-2018: Modify script to automatically make plots directory
#                     and save.
#
#       (Old stuff: I don't remember the dates (08/2017-08/2018))
#       * Change functions to adapt Julia 1.0 changes
#       * Move the cumulative slip calculation outside the time loop 
#       * Add scripts to compute the earthquake magnitude     
#       * Add plots script for various plotting functions
#       * Implemented elastic halfspace
#       * Setup for a shallow fault zone
#       * Organize everything into structs and functions
#       * Interpolation for initial stress and friction in halfspace
#       * Add separate files for parameters, initial conditions, functions
###############################################################################

# Output results
mutable struct results
    Stress::Array{Float64,2}
    SlipVel::Array{Float64,2}
    Slip::Array{Float64,2}
    time_::Array{Float64}
end


function main(P::parameters, S::input_variables)

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
    d = SharedArray{Float64}(S.nglob)
    v = SharedArray{Float64}(S.nglob)
    v .= 0.5e-3
    a = SharedArray{Float64}(S.nglob)
    
    #.....................................
    # Stresses and time related variables
    #.....................................
    tau = SharedArray{Float64}(P.FltNglob)
    FaultC = SharedArray{Float64}(P.FltNglob)
    Vf = SharedArray{Float64}(P.FltNglob)
    Vf1 = SharedArray{Float64}(P.FltNglob)
    Vf2 = SharedArray{Float64}(P.FltNglob)
    Vf0 = SharedArray{Float64}(length(S.iFlt))
    FltVfree = SharedArray{Float64}(length(S.iFlt))
    psi = SharedArray{Float64}(P.FltNglob)
    psi0 = SharedArray{Float64}(P.FltNglob)
    psi1 = SharedArray{Float64}(P.FltNglob)
    psi2 = SharedArray{Float64}(P.FltNglob)
    tau1 = SharedArray{Float64}(P.FltNglob)
    tau2 = SharedArray{Float64}(P.FltNglob)
    tau3 = SharedArray{Float64}(P.FltNglob)


    # Initial state variable
    psi = S.tauo./(S.Seff.*S.ccb) - P.fo./S.ccb - (S.cca./S.ccb).*log.(2*v[S.iFlt]./P.Vo)
    psi0 .= psi[:]

    isolver::Int = 1
    
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

    F = SharedArray{Float64}(S.nglob)
    dPre = SharedArray{Float64}(S.nglob)
    vPre = SharedArray{Float64}(S.nglob)
    dd = SharedArray{Float64}(S.nglob)
    dacum = SharedArray{Float64}(S.nglob)
    dnew = SharedArray{Float64}(length(S.FltNI))

    # Preallocate variables with unknown size
    output = results(zeros(P.FltNglob, 1000000), zeros(P.FltNglob, 1000000), 
                         zeros(P.FltNglob, 1000000), zeros(1000000))
    # Iterators
    idelevne = 3
    tevneb = 0
    tevne = 0

    v = v[:] .- 0.5*P.Vpl
    Vf = 2*v[S.iFlt]
    iFBC::Array{Int64} = findall(abs.(S.FltX) .> 24e3)
    NFBC::Int64 = length(iFBC)
    Vf[iFBC] .= 0


    v[S.FltIglobBC] .= 0

    #....................
    # Start of time loop
    #....................
    it = 0
    t = 0

    while t < P.Total_time
        it = it + 1
        t = t + dt

        output.time_[it] = t 

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
                              S.H, S.Ht, S.iglob, S.nglob, S.W)
                
                # update displacement on the medium
                d[S.FltNI] .= dnew

                # make d = F on the fault
                d[S.iFlt] .= F[S.iFlt]

                # Compute on-fault stress
                a .= 0

                # Compute forcing (acceleration) for each element
                #println("\nElement Computation:")
                a = element_computation(P, S.iglob, d, S.H, S.Ht, S.W, a)

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
            a = element_computation2(P, S.iglob, d, S.H, S.Ht, S.W, a)
            a[S.FltIglobBC] .= 0

            # Absorbing boundaries
            a[S.iBcL] .= a[S.iBcL] .- S.BcLC.*v[S.iBcL]
            a[S.iBcT] .= a[S.iBcT] .- S.BcTC.*v[S.iBcT]

            ###### Fault Boundary Condition: Rate and State #############
            FltVfree .= 2*v[S.iFlt] .+ 2*half_dt*a[S.iFlt]./S.M[S.iFlt]
            Vf .= 2*vPre[S.iFlt] .+ P.Vpl


            # Sliprate and NR search
            psi1, Vf1, tau1, psi2, Vf2, tau2 = FBC(P, S, NFBC, psi1, Vf1, 
                                    tau1, psi2, Vf2, tau2, psi, Vf, FltVfree, dt)
            

            tau .= tau2 .- S.tauo
            tau[iFBC] .= 0
            psi .= psi2
            #KD = a[:]
            a[S.iFlt] .= a[S.iFlt] .- S.FltB.*tau
            ########## End of fault boundary condition ############## 


            # Solve for a_new
            a .= a./S.M
            
            # Correction
            v .= v .+ half_dt*a

            v[S.FltIglobBC] .= 0
            a[S.FltIglobBC] .= 0

            #### Line 861: Omitting P_Ma
            
        end # of isolver if loop
        
        Vfmax = 2*maximum(v[S.iFlt]) .+ P.Vpl


        #----
        # Output variables at different depths for every timestep
        # Omitted the part of code from line 871 - 890, because I 
        # want to output only certain variables each timestep
        # Doing it in separate script
        #----


        #-----
        # Output stress and slip before and after events
        # Doing it in separate script
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
        output.Stress[:,it] = (tau + S.tauo)./1e6
        output.SlipVel[:,it] = 2*v[S.iFlt] .+ P.Vpl
        output.Slip[:,it] = 2*d[S.iFlt] .+ P.Vpl*t
        
        # Compute next timestep dt
        dt = dtevol(P, dt , dtmin, S.XiLf, P.FltNglob, NFBC, output.SlipVel[:,it], isolver)

    end # end of time loop
    
    # Remove zeros from preallocated vectors
    output.time_ = output.time_[1:it]
    output.Stress = output.Stress[:,1:it]
    output.SlipVel = output.SlipVel[:,1:it]
    output.Slip = output.Slip[:,1:it]

    return output #results(Stress, SlipVel, Slip, time_)

end
