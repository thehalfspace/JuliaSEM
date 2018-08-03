# Calculate XiLf used in computing the timestep

function XiLfFunc(s::space_parameters, t::time_parameters, 
                  eq::earthquake_parameters, muMax, cca, ccb, Seff)

    hcell = s.LX/(s.FltNglob-1)
    Ximax = 0.5
    Xithf = 1

    Xith = Array{Float64}(s.FltNglob)
    XiLf = Array{Float64}(s.FltNglob)

    for j = 1:s.FltNglob

        # Compute time restricting parameters
        expr1 = -(cca[j] - ccb[j])/cca[j]
        expr2 = t.gamma_*muMax/hcell*eq.xLf[j]/(cca[j]*Seff[j])
        ro = expr2 - expr1

        if (0.25*ro*ro - expr2) >= 0
            Xith[j] = 1/ro
        else
            Xith[j] = 1 - expr1/expr2
        end

        # For each node, compute slip that node cannot exceed in one timestep
        if Xithf*Xith[j] > Ximax
            XiLf[j] = Ximax*eq.xLf[j]
        else
            XiLf[j] = Xithf*Xith[j]*eq.xLf[j]
        end
    end
       

    return XiLf
end
