###################################
#      SOLVER TIME LOOP           #
# Explicit Newmark Scheme with    #
# alpha=1, beta=0, gamma=0.5      #
###################################

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

diagKnew - Kdiag[FltNI]


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

        end

    end

end
