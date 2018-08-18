###################################################
#	ASSEMBLE THE MASS AND THE STIFFNESS MATRICES
###################################################


    
function assemble(P::parameters, iglob, M, W)


    xgll, wgll, H = GetGLL(P.NGLL)
    wgll2 = wgll*wgll';

    vso = zeros(P.NGLL, P.NGLL)
    vs = zeros(P.NGLL-1, P.NGLL)
    dx = zeros(P.NGLL-1, P.NGLL)
    muMax = 0
    dt = Inf

    for ey = 1:P.NelY
        for ex = 1:P.NelX

            eo = (ey-1)*P.NelX + ex
            ig = iglob[:,:,eo]

            # Properties of heterogeneous medium
            if ex*P.dxe >= P.ThickX && ey*P.dye <= P.ThickY
                P.rho[:,:] .= P.rho2
                P.mu[:,:] .= P.rho2*P.vs2^2
            else
                P.rho[:,:] .= P.rho1
                P.mu[:,:] .= P.rho1*P.vs1^2
            end

            if muMax < maximum(maximum(P.mu))
                muMax = maximum(maximum(P.mu))
            end

            # Diagonal Mass Matrix
            M[ig] .+= wgll2.*P.rho*P.jac

            # Local contributions to the stiffness matrix
            W[:,:,eo] .= wgll2.*P.mu;
            
            # Set timestep
            vso .= sqrt.(P.mu./P.rho)
            
            if P.dxe<P.dye
                vs .= max.(vso[1:P.NGLL-1,:], vso[2:P.NGLL,:])
                dx .= repeat( diff(xgll)*0.5*P.dxe, 1, P.NGLL)
            else
                vs .= max.(vso[:,1:P.NGLL-1], vso[:,2:P.NGLL])'
                dx .= repeat( diff(xgll)*0.5*P.dye, 1, P.NGLL)
            end
            
            dtloc = dx./vs
            dt = minimum( push!(dtloc[1:end], dt) )

        end
    end

    return M, W, dt, muMax
end
