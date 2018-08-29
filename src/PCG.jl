################################################
#                                              
#   SOLVE FOR DISPLACEMENT USING PRECONDITIONED 
#           CONJUGATE GRADIENT METHOD          
#                                              
################################################
#using JLD2

function PCG(P::parameters, diagKnew, dnew, F, iFlt,
             FltNI, H, Ht, iglob, nglob, W)
    
    a_local = SharedArray{Float64}(nglob)
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

            #filename = string(dir, "/data", name, "pcgfail.jld2")
            #@save filename dnew rnew Fnew
            @error("PCG did not converge")
            return
        end
    end

    return dnew
end


# Sub function to be used inside PCG
function element_computation(P::parameters, iglob, F_local, H, Ht, W, a_local)
    
    @sync @distributed for eo = 1:P.Nel

        ig = iglob[:,:,eo]

        locall = local_calc(P, F_local[ig], H, Ht, W[:,:,eo])
        
        #  println(locall)

        # Assemble into global vector
        a_local[ig] .+= locall

    end

    return a_local

end

@everywhere function local_calc(P, locall, H, Ht, Wlocal)

    # Gradients wrt local variables
    d_xi = Ht*locall
    d_eta = locall*H

    # Element contribution to the internal forces
    locall = P.coefint1*H*(Wlocal.*d_xi) + 
            P.coefint2*(Wlocal.*d_eta)*Ht
    return locall 
end

# Compute the displacement/forcing for each element
function element_computation2(P::parameters, iglob, F_local, H, Ht, W, a_local)
    
   @sync @distributed for eo = 1:P.Nel

        # Switch to local element representation
        ig = iglob[:,:,eo]

        locall = local_calc(P, F_local[ig], H, Ht, W[:,:,eo])

        # Assemble into global vector
        a_local[ig] .-= locall

    end

    return a_local

end
