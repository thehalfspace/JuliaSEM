###################################
#      SOLVER TIME LOOP           #
# Explicit Newmark Scheme with    #
# alpha=1, beta=0, gamma=0.5      #
###################################

include("PCG.jl")
include("otherFunctions.jl")

# Initially 1 = quasistatic; 2 = dynamic
isolver = 1

# Compute the diagonal of K
Kdiag = zeros(nglob,1)
Klocdiag = zeros(NGLL, NGLL)
for e = 1:Nel
    ig = iglob[:,:,e]
    wloc = W[:,:,e]
    Klocdiag[:,:] = 0

    for k =  1:NGLL
        for j = 1:NGLL
            Klocdiag[k,j] = Klocdiag[k,j] + 
                            sum(coefint1*H[k,:]'.*(wloc[:,j].*Ht[:,k])
                            + coefint2*(wloc[k,:].*H[j,:])'.*Ht[:,j])
        end
    end

    Kdiag[ig] = Kdiag[ig] + Klocdiag[:,:]
end

diagKnew = Kdiag[FltNI]


v[:] = v[:] - 0.5*Vpl
Vf = 2*v[iFlt]
iFBC = find(abs.(FltX) .> 24e3/distN)
NFBC = length(iFBC)
Vf[iFBC] = 0

FltIglobBC = iFlt[iFBC] # Fault boundary

v[FltIglobBC] = 0



# Output sliprate at the start of every cycle
flag = 0
event_iter1 = 1
event_iter2 = 1


# Preallocate variables with unknown size
time_ = [0]

#...........................
# START OF THE TIME LOOP
#...........................

while t < 10 #Total_time
    it = it + 1
    t = t + 1#dt

    time_ = push!(time_,t) 


    if isolver == 1

        vPre = v
        dPre = d

        Vf0 = 2*v[iFlt] + Vpl
        Vf = Vf0

        for p1 = 1:2
            
            # Compute the forcing term
            F[:] = 0
            F[iFlt] = dPre[iFlt] + v[iFlt]*dt

            # Assign previous displacement field as initial guess
            dnew = d[FltNI]

            # Solve d = K^-1F by PCG
            dnew = PCG(coefint1, coefint2, diagKnew, dnew, F, iFlt, FltNI,
                          H, Ht, iglob, Nel, nglob, W, a)
            
            # update displacement on the medium
            d[FltNI] = dnew

            # make d = F on the fault
            d[iFlt] = F[iFlt]

            # Compute on-fault stress
            a[:] = 0
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

            a[FltIglobBC] = 0
            tau1 = -a[iFlt]./FltB
            
            # Compute slip-rates on-fault
            for jF = 1:FaultNglob-NFBC 

                j = jF - 1 + NFBC
                psi1[j] = IDS(psi[j], psi1[j], dt, Vo[j], xLf[j], Vf[j])

                tauAB[j] = tau1[j] + tauo[j]
                fa = tauAB[j]/(Seff[j]*cca[j])
                help = -(fo[j] + ccb[j]*psi1[j])/cca[j]
                help1 = exp(help + fa)
                help2 = exp(help - fa)
                Vf1[j] = Vo[j]*(help1 - help2) 
            end
        end

        psi = psi1
        tau = tau1
        tau[iFBC] = 0
        Vf1[iFBC] = Vpl

        v[iFlt] = 0.5*(Vf1 - Vpl)
        v[FltNI] = (d[FltNI] - dPre[FltNI])/dt

        RHS = a
        RHS[iFlt] = RHS[iFlt] - FltB.*tau
        RMS = sqrt(sum(RHS.^2)/length(RHS))./maximum(abs.(RHS))
        
        # Line 731: P_MA: Omitted
        a[:] = 0
        d[FltIglobBC] = 0
        v[FltIglobBC] = 0

        
        # If isolver != 1, or max slip rate is < 10^-2 m/s
    else

        dPre = d
        vPre = v

        # Update
        d = d + dt*v + (half_dt^2)*a

        # Prediction
        v = v + half_dt*a
        a[:] = 0


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
        
            # Element contribution to the internal forces
            locall = coefint1*H*(W[:,:,eo].*d_xi) + 
                    coefint2*(W[:,:,eo].*d_eta)*Ht

            # Assemble into global vector
            a[ig] = a[ig] + locall
        end

        a[FltIglobBC] = 0

        # Absorbing boundaries
        a[iBcL] = a[iBcL] - BcLC.*v[iBcL]
        a[iBcT] = a[iBcT] - BcTC.*v[iBcT]

        ###### Fault Boundary Condition: Rate and State #############
        FltVfree = 2*v[iFlt] + 2*half_dt*a[iFlt]./M[iFlt]
        Vf = 2*vPre[iFlt] + Vpl

        for jF = 1:FaultNglob-NFBC

            j = jF - 1 + NFBC
            psi1[j] = IDS(psi[j], psi1[j], dt, Vo[j], xLf[j], Vf[j])

            Vf1[j], tau1[j] = NRsearch(fo[j], Vo[j], cca[j], ccb[j],Seff[j],
                                      tauNR[j], tauo[j], psi1[j], FltZ[j], FltVfree[j])
        
            if Vf[j] > 1e10 || isnan(Vf[j]) == 1 || isnan(tau1[j]) == 1
                error("NR SEARCH FAILED!")
                return
            end
             
            if IDstate == 1
                psi2[j] = psi[j] + 0.5*dt*( (Vo[j]/xLf[j])*exp(-psi[j]) 
                                            - abs(Vf[j])/xLf[j] 
                                            + (Vo[j]/xLf[j])*exp(-psi1[j]) 
                                            - abs(Vf1[j])/xLf[j] )

            elseif IDstate == 2
                VdtL = 0.5*abs(Vf1[j] + Vf[j])*dt/xLf[j]
                if VdtL < 1e-6
                    psi2[j] = log( exp(psi[j]-VdtL) + Vo[j]*dt/xLf[j] -
                                  0.5*Vo[j]*0.5*abs(Vf1[j] + Vf[j])*dt*dt/(xLf[j]^2))
                else
                    psi2[j] = log(exp(psi[j]-VdtL) + 
                                  (Vo[j]/(0.5*abs(Vf[j] + Vf1[j])))*(1-exp(-VdtL)))
                end

            elseif IDstate == 3
                psi2[j] = exp(-0.5*abs(Vf[j] + Vf1[j])*dt/xLf[j]) * 
                                log(0.5*abs(Vf[j] + Vf1[j])/Vo[j]) + 
                                exp(-0.5*abs(Vf[j] + Vf1[j])*dt/xLf[j])*psi[j] 
                                + log(Vo[j]/(-0.5*abs(Vf[j] + Vf1[j])) )

                if ~any(imag(psi1)) == 0
                    return
                end
            end

            # NRsearch 2nd loop
            Vf2[j], tau2[j] = NRsearch(fo[j], Vo[j], cca[j], ccb[j],Seff[j],
                                      tau1[j], tauo[j], psi2[j], FltZ[j], FltVfree[j])

        end
        
        tau = tau2 - tauo
        tau[iFBC] = 0
        psi = psi2
        KD = a
        a[iFlt] = a[iFlt] - FltB.*tau
        ########## End of fault boundary condition ############## 


        RHS = a

        # Solve for a_new
        a = a./M
        
        # Correction
        v = v + half_dt*a

        v[FltIglobBC] = 0
        a[FltIglobBC] = 0

        #### Line 861: Omitting P_Ma
        
        LHS = M.*a
        RMS = sqrt(sum.((RHS - LHS).^2)/length(RHS))./max(abs(RHS))

    end # of isolver if loop
    
    Vfmax = 2*maximum(v[iFlt]) + Vpl

    #----
    # Output variables at different depths for every timestep
    # Omitted the part of code from line 871 - 890, because I 
    # want to output only certain variables each timestep
    #----

    # Output stress, slip, sliprate on fault every certain interval
    # line 891


end


