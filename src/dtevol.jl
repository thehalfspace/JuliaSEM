############################################
#   Compute the timestep for next iteration
############################################

function dtevol(P::parameters, dt, dtmin, XiLf, FaultNglob, NFBC, Vf, isolver)

    if isolver == 1

        # initial value of dt
        dtnx = P.dtmax

        # Adjust the timestep according to cell velocities and slip
        for i = NFBC:FaultNglob - 1

            if abs(Vf[i])*P.dtmax > XiLf[i]
                dtcell = XiLf[i]/abs(Vf[i])

                if dtcell < dtnx
                    dtnx = dtcell
                end
            end
        end

        if dtmin > dtnx
            dtnx = dtmin
        end

        if dtnx > P.dtincf*dt
            dtnx = P.dtincf*dt
        end

        dt = dtnx

    elseif isolver == 2
        
        dt = dtmin
    end

    return dt

end
