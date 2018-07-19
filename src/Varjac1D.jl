##########################
# Variable mesh jacobian
##########################

function Varjac1D(dxe, dye)

    dx_dxi = 0.5.*dxe
    dy_deta = 0.5.*deta

    return dx_dxi, dy_deta, jac, coefint1, coefint2
end
