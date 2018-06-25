###################################
#      SOLVER TIME LOOP           #
# Explicit Newmark Scheme with    #
# alpha=1, beta=0, gamma=0.5      #
###################################

include("PCG.jl")

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

while t < 10#Total_time
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
                          H, Ht, iglob, Nel, nglob, W)
            
            # update displacement on the medium
            d[FltNI] = dnew

            # make d = F on the fault
            d[iFlt] = f[iFlt]

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
        end

    end

end
