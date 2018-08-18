################################################
#                                              
#   SOLVE FOR DIPLACEMENT USING PRECONDITIONED 
#           CONJUGATE GRADIENT METHOD          
#                                              
################################################

using JLD2

function PCG(P::parameters, diagKnew, dnew, F, iFlt,
             FltNI, H, Ht, iglob, nglob, W)
    
    #global a 
    a_local = zeros(nglob)
    dd_local = zeros(nglob)
    p_local = zeros(nglob)
    
    a_local = element_computation(P, iglob, F, H, Ht, W, a_local)
    Fnew = -a_local[FltNI]

    dd_local[FltNI] .= dnew
    dd_local[iFlt] .= 0

    a_local[:] .= 0
    
    a_local = element_computation(P, iglob, dd_local, H, Ht, W, a_local)
    anew = a_local[FltNI]

    # Initial residue
    rnew = Fnew - anew
    znew = rnew./diagKnew
    pnew = znew
    p_local[:] .= 0
    p_local[FltNI] = pnew

    @inbounds for n = 1:4000
        anew[:] .= 0
        a_local[:] .= 0
        
        a_local = element_computation(P, iglob, p_local, H, Ht, W, a_local)

        anew = a_local[FltNI]

        alpha_ = znew'*rnew/(pnew'*anew)
        dnew  .+= alpha_*pnew
        rold = rnew
        zold = znew
        rnew = rold - alpha_*anew
        znew = rnew./diagKnew
        beta_ = znew'*rnew/(zold'*rold)
        pnew = znew + beta_*pnew
        p_local[:] .= 0
        p_local[FltNI] = pnew

        if norm(rnew)/norm(Fnew) < 1e-5
            break;
        end

        if n == 4000 || norm(rnew)/norm(Fnew) > 1e10
            print(norm(rnew)/norm(Fnew))
            println("n = ", n)

            filename = string(dir, "/data", name, "pcgfail.jld")
            @save filename dnew, rnew, Fnew
            @error("PCG did not converge")
            return
        end
    end

    return dnew
end


# Sub function to be used inside PCG
function element_computation(P::parameters, iglob, F_local, H, Ht, W, a_local)
    
    @inbounds for eo = 1:P.Nel

        # Switch to local element representation
        ig = iglob[:,:,eo]
        locall = F_local[ig]

        # Gradients wrt local variables
        d_xi = Ht*locall
        d_eta = locall*H

        # Element contribution to the internal forces
        locall = P.coefint1*H*(W[:,:,eo].*d_xi) + 
                 P.coefint2*(W[:,:,eo].*d_eta)*Ht

        # Assemble into global vector
        a_local[ig] += locall

    end

    return a_local

end

# Compute the displacement/forcing for each element
function element_computation2(P::parameters, iglob, F_local, H, Ht, W, a_local)
    
    @inbounds for eo = 1:P.Nel

        # Switch to local element representation
        ig = iglob[:,:,eo]
        locall = F_local[ig]

        # Gradients wrt local variables
        d_xi = Ht*locall
        d_eta = locall*H

        # Element contribution to the internal forces
        locall = P.coefint1*H*(W[:,:,eo].*d_xi) + 
                 P.coefint2*(W[:,:,eo].*d_eta)*Ht

        # Assemble into global vector
        a_local[ig] -= locall

    end

    return a_local

end
