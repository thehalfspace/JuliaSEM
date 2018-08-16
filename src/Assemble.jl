###################################################
#	ASSEMBLE THE MASS AND THE STIFFNESS MATRICES
###################################################


    
function assemble(s::space_parameters, m::medium_properties, iglob, M, W)


    xgll, wgll, H = GetGLL(s.NGLL)
    wgll2 = wgll*wgll';

    vso = zeros(s.NGLL, s.NGLL)
    vs = zeros(s.NGLL-1, s.NGLL)
    dx = zeros(s.NGLL-1, s.NGLL)
    muMax = 0
    dt = Inf

    for ey = 1:s.NelY
        for ex = 1:s.NelX

            eo = (ey-1)*s.NelX + ex
            ig = iglob[:,:,eo]

            # Properties of heterogeneous medium
            if ex*s.dxe >= m.ThickX && ey*s.dye <= m.ThickY
                m.rho[:,:] .= m.rho2
                m.mu[:,:] .= m.rho2*m.vs2^2
            else
                m.rho[:,:] .= m.rho1
                m.mu[:,:] .= m.rho1*m.vs1^2
            end

            if muMax < maximum(maximum(m.mu))
                muMax = maximum(maximum(m.mu))
            end

            # Diagonal Mass Matrix
            M[ig] .+= wgll2.*m.rho*s.jac

            # Local contributions to the stiffness matrix
            W[:,:,eo] .= wgll2.*m.mu;
            
            # Set timestep
            vso .= sqrt.(m.mu./m.rho)
            
            if s.dxe<s.dye
                vs .= max.(vso[1:s.NGLL-1,:], vso[2:s.NGLL,:])
                dx .= repeat( diff(xgll)*0.5*s.dxe, 1, s.NGLL)
            else
                vs .= max.(vso[:,1:s.NGLL-1], vso[:,2:s.NGLL])'
                dx .= repeat( diff(xgll)*0.5*s.dye, 1, s.NGLL)
            end
            
            dtloc = dx./vs
            dt = minimum( push!(dtloc[1:end], dt) )

        end
    end

    return M, W, dt, muMax
end
