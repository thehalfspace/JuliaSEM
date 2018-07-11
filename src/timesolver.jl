###################################
#      SOLVER TIME LOOP           #
# Explicit Newmark Scheme with    #
# alpha=1, beta=0, gamma=0.5      #
###################################

include("PCG.jl")
include("dtevol.jl")
include("NRsearch.jl")
#include("otherFunctions.jl")

temp = [0.0]
Vf0 = zeros(length(iFlt))
FltVfree = zeros(length(iFlt))

#...........................
# START OF THE TIME LOOP
#...........................

while it < 100
    it = it + 1
    t = t + dt

    time_[it] = t 


    if isolver == 1

        vPre .= v[:]
        dPre .= d[:]

        Vf0 .= 2*v[iFlt] + Vpl
        Vf  .= Vf0[:]

        for p1 = 1:2
            
            # Compute the forcing term
            F[:] .= 0
            F[iFlt] .= dPre[iFlt] .+ v[iFlt]*dt

            # Assign previous displacement field as initial guess
            dnew .= d[FltNI]

            # Solve d = K^-1F by PCG
            dnew = PCG(coefint1, coefint2, diagKnew, dnew, F, iFlt, FltNI,
                          H, Ht, iglob, Nel, nglob, W)
            
            # update displacement on the medium
            d[FltNI] .= dnew

            # make d = F on the fault
            d[iFlt] .= F[iFlt]

            # Compute on-fault stress
            a[:] .= 0
            for eo = 1:Nel
                ig = iglob[:,:,eo]
                locall = d[ig]

                # Gradients wrt local variables
                d_xi = Ht*locall
                d_eta = locall*H

                # Element contribution
                wloc = W[:,:,eo]
                d_xi = wloc.*d_xi
                d_xi = H*d_xi
                d_eta = wloc.*d_eta
                d_eta = d_eta*Ht
                locall = coefint1*d_xi + coefint2*d_eta

                # Assemble into global vector
                a[ig] = a[ig] + locall
            end

            a[FltIglobBC] .= 0
            tau1 .= -a[iFlt]./FltB
            
            # Compute slip-rates on-fault
            for jF = 1:FaultNglob-NFBC 

                j = jF - 1 + NFBC
                psi1[j] = IDS(psi[j], dt, Vo[j], xLf[j], Vf[j], 1e-6, IDstate)

                tauAB[j] = tau1[j] + tauo[j]
                fa = tauAB[j]/(Seff[j]*cca[j])
                help = -(fo[j] + ccb[j]*psi1[j])/cca[j]
                help1 = exp(help + fa)
                help2 = exp(help - fa)
                Vf1[j] = Vo[j]*(help1 - help2) 
            end
            
            Vf1[iFBC] .= Vpl
            Vf .= (Vf0 + Vf1)/2
            v[iFlt] .= 0.5*(Vf - Vpl)

        end

        psi .= psi1[:]
        tau .= tau1[:]
        tau[iFBC] .= 0
        Vf1[iFBC] .= Vpl

        v[iFlt] .= 0.5*(Vf1 - Vpl)
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
        for eo = 1:Nel
            
            # Switch to local element representation
            ig = iglob[:,:,eo]
            isETA = eo<=Nel_ETA
            if isETA
                locall = d[ig] + eta.*v[ig]
            else
                locall = d[ig]
            end
            # Gradients wrt local variables
            d_xi = Ht*locall
            d_eta = locall*H
        
            # Element contribution
            wloc = W[:,:,eo]
            d_xi = wloc.*d_xi
            d_xi = H*d_xi
            d_eta = wloc.*d_eta
            d_eta = d_eta*Ht
            locall = coefint1*d_xi + coefint2*d_eta

            # Assemble into global vector
            a[ig] = a[ig] - locall
        end

        a[FltIglobBC] .= 0

        # Absorbing boundaries
        a[iBcL] .= a[iBcL] .- BcLC.*v[iBcL]
        a[iBcT] .= a[iBcT] .- BcTC.*v[iBcT]

        ###### Fault Boundary Condition: Rate and State #############
        FltVfree .= 2*v[iFlt] .+ 2*half_dt*a[iFlt]./M[iFlt]
        Vf .= 2*vPre[iFlt] .+ Vpl

        for jF = 1:FaultNglob-NFBC

            j = jF - 1 + NFBC
            psi1[j] = IDS(psi[j], dt, Vo[j], xLf[j], Vf[j], 1e-5, IDstate)

            Vf1[j], tau1[j] = NRsearch(fo[j], Vo[j], cca[j], ccb[j],Seff[j],
                                      tauNR[j], tauo[j], psi1[j], FltZ[j], FltVfree[j])
        
            if Vf[j] > 1e10 || isnan(Vf[j]) == 1 || isnan(tau1[j]) == 1
                error("NR SEARCH FAILED!")
                return
            end
            
            psi2[j] = IDS2(psi[j], psi1[j], dt, Vo[j], xLf[j], Vf[j], Vf1[j], IDstate)
            
            # NRsearch 2nd loop
            Vf2[j], tau2[j] = NRsearch(fo[j], Vo[j], cca[j], ccb[j],Seff[j],
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
    
    Vfmax = 2*maximum(v[iFlt]) + Vpl

    #----
    # Output variables at different depths for every timestep
    # Omitted the part of code from line 871 - 890, because I 
    # want to output only certain variables each timestep
    #----

    # Output stress, slip, sliprate on fault every certain interval
    if t>tvsx
        ntvsx = ntvsx + 1
        
        delf5yr[:,ntvsx] = 2*d[iFlt] + Vpl*t
        Vf5yr[:,ntvsx] = 2*v[iFlt] + Vpl
        Tau5yr[:,ntvsx] = (tau + tauo)./1e6
        
        tvsx = tvsx + tvsxinc
    end
    
    if Vfmax > Vevne 
        if idelevne == 0
            nevne = nevne + 1
            idelevne = 2
            tevneb = t
            tevne = tevneinc

            delfsec[:,nevne] = 2*d[iFlt] + Vpl*t
            Vfsec[:,nevne] = 2*v[iFlt] + Vpl
            Tausec[:,nevne] = (tau + tauo)./1e6
        end
        if idelevne == 1 && (t - tevneb) > tevne
            nevne = nevne + 1
            
            delfsec[:,nevne] = 2*d[iFlt] + Vpl*t
            Vfsec[:,nevne] = 2*v[iFlt] + Vpl
            Tausec[:,nevne] = (tau + tauo)./1e6

            tevne = tevne + tevneinc
        end

    else
        global idelevne = 0
    end

    #-----
    # Output stress and slip before and after events
    # Omitting lines 920-934
    #-----

    # Output timestep info on screen
    if mod(it,40) == 0
        @printf("\nTime (yr) = %1.5g", t/yr2sec)
    end
    
    # Determine quasi-static or dynamic regime based on max-slip velocity
    if isolver == 1 && Vfmax < 5e-3 || isolver == 2 && Vfmax < 2e-3
        isolver = 1
    else
        isolver = 2
    end

    # Some variables for each timestep
    Stress[:,it] = (tau + tauo)./1e6
    SlipVel[:,it] = 2*v[iFlt] + Vpl
    Slip[:,it] = 2*d[iFlt] + Vpl*t

    temp = push!(temp,Vf[120])
    
    # Compute next timestep dt
    dt = dtevol(dt, dtmax, dtmin, dtincf, XiLf, FaultNglob, NFBC, SlipVel[:,it], isolver)


end # End of time loop

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

println("\nSimulation Complete")
