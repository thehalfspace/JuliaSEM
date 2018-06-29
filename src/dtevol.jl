############################################
#   Compute the timestep for next iteration
############################################

function dtevol(dt, dtmax, dtmin, dtincf, XiLf, FaultNglob, NFBC, Vf, isolver)

    if isolver == 1

        # initial value of dt
        dtnx = dtmax

        # Adjust the timestep according to cell velocities and slip
        for iF = 1:FaultNglob - NFBC
            i = iF - 1 + NFBC

            if abs(Vf[i])*dtmax > XiLf[i]
                dtcell = XiLf[i]/abs(Vf[i])

                if dtcell < dtnx
                    dtnx = dtcell
                end
            end
        end

        if dtmin > dtnx
            dtnx = dtmin
        end

        if dtnx > dtincf*dt
            dtnx = dtincf*dt
        end

        dt = dtnx

    elseif isolver == 2
        
        dt = dtmin
    end

    return dt

end
