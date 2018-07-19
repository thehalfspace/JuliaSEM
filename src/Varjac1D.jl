##########################
# Variable mesh jacobian
##########################

#function Varjac1D(dxe, dye)

#    dx_dxi = 0.5.*dxe
#    dy_deta = 0.5.*deta

#    return dx_dxi, dy_deta, jac, coefint1, coefint2
#end


function coefint(coefint1, coefint2, dxe, dye)
   
    e = 0
    for ey=1:length(dye)
        for ex=1:length(dxe)
            
            e = e+1
            jac = 0.25*dxe[ex]*dye[ey]
    
            coefint1[e] = jac/(0.5*dxe[ex])^2
            coefint2[e] = jac/(0.5*dye[ey])^2
        end
    end

    return coefint1, coefint2
end
