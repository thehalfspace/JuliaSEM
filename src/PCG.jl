################################################
#                                              #
#   SOLVE FOR DIPLACEMENT USING PRECONDITIONED #
#           CONJUGATE GRADIENT METHOD          #
#                                              #
################################################


function PCG(coefint1, coefint2, diagKnew, dnew, F, iFlt,
             FltNI, H, Ht, iglob, Nel, nglob, W)
    
    
    a = zeros(nglob)
    dd = zeros(nglob)
    p = zeros(nglob)
    
    a = element_computation(Nel, iglob, F, H, Ht, coefint1, coefint2, W)
    Fnew = -a[FltNI]

    dd[FltNI] = dnew
    dd[iFlt] = 0

    a[:] = 0
    
    a = element_computation(Nel, iglob, dd, H, Ht, coefint1, coefint2, W)
    anew = a[FltNI]

    # Initial residue
    rnew = Fnew - anew
    znew = rnew./diagKnew
    pnew = znew
    p[:] = 0
    p[FltNI] = pnew

    for n = 1:4000
        anew[:] = 0
        a[:] = 0
        
        a = element_computation(Nel, iglob, p, H, Ht, coefint1, coefint2, W)

        anew = a[FltNI]
        alpha_ = znew'*rnew/(pnew'*anew)
        dnew = dnew + alpha_*pnew
        rold = rnew
        zold = znew
        rnew = rold - alpha_*anew
        znew = rnew./diagKnew
        beta_ = znew'*rnew/(zold'*rold)
        pnew = znew + beta_*pnew
        p[:] = 0
        p[FltNI] = pnew

        if norm(rnew)/norm(Fnew) < 1e-5
            break;
        end

        if n == 4000 || norm(rnew)/norm(Fnew) > 1e10
            println("n = ", n)
            error("PCG did not converge")
            return
        end
    end

    return dnew
end


# Sub function to be used inside PCG
function element_computation(Nel, iglob, F, H, Ht, coefint1, coefint2, W)
    
    for eo = 1:Nel

        # Switch to local element representation
        ig = iglob[:,:,eo]
        locall = F[ig]

        # Gradients wrt local variables
        d_xi = Ht*locall
        d_eta = locall*H

        # Element contribution to the internal forces
        locall = coefint1*H*(W[:,:,eo].*d_xi) + 
                 coefint2*(W[:,:,eo].*d_eta)*Ht

        # Assemble into global vector
        a[ig] = a[ig] + locall

    end

    return a

end
