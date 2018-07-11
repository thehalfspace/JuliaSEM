###################################################
#	ASSEMBLE THE MASS AND THE STIFFNESS MATRICES
###################################################


    
rho = zeros(NGLL, NGLL)
mu = zeros(NGLL, NGLL)
vso = zeros(NGLL, NGLL)
vs = zeros(NGLL-1, NGLL)
dx = zeros(NGLL-1, NGLL)
muMax = 0

for ey = 1:NelY
    for ex = 1:NelX

        eo = (ey-1)*NelX + ex
        ig = iglob[:,:,eo]

        # Properties of heterogeneous medium
        if ex*dxe >= ThickX && ey*dye <= ThickY
            rho[:,:] .= rho2
            mu[:,:] .= rho2*vs2^2
        else
            rho[:,:] .= rho1
            mu[:,:] .= rho1*vs1^2
        end

        if muMax < maximum(maximum(mu))
            muMax = maximum(maximum(mu))
        end

        # Diagonal Mass Matrix
        M[ig] .= M[ig] .+ wgll2.*rho*jac

        # Local contributions to the stiffness matrix
        W[:,:,eo] .= wgll2.*mu;
            
        # Set timestep
        vso .= sqrt.(mu./rho)
            
        if dxe<dye
            vs .= max.(vso[1:NGLL-1,:], vso[2:NGLL,:])
            dx .= repmat( diff(xgll)*0.5*dxe, 1, NGLL)
        else
            vs .= max.(vso[:,1:NGLL-1], vso[:,2:NGLL])'
            dx .= repmat( diff(xgll)*0.5*dye, 1, NGLL)
        end
            
        dtloc = dx./vs
        dt = minimum( push!(dtloc[1:end], dt) )


    end
end


#=
# With single for loop

for c = 1:Nel

    ig = iglob[:,:,c]

    # Properties of heterogeneous medium
    rho[:,:] = 1
    mu[:,:] = 1

    # Diagonal Mass Matrix
    M[ig] = M[ig] + wgll2.*rho*jac

    # Local contributions to the stiffness matrix
    W[:,:,c] = wgll2.*mu;
        
    # Set timestep
    vs = sqrt.(mu./rho)

    if dxe<dye
        vs = max.(vs[1:NGLL-1,:], vs[2:NGLL,:])
        dx = repmat( diff(xgll)*0.5*dxe, 1, NGLL)
    else
        vs = max.(vs[:,1:NGLL-1], vs[:,2:NGLL])
        dx = repmat( diff(xgll)*0.5*dye, 1, NGLL)'
    end
        
    dtloc = dx./vs
    dt = minimum( push!(dtloc[1:end], dt) )    

end
=#;
