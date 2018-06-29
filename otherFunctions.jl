#########################
# SOME OTHER FUNCTIONS
#########################

function IDS(psi, psi1, dt, Vo, xLf, Vf, IDstate = 2)
    #= compute slip-rates on fault based on different
       formulations =#

    if IDstate == 1
        psi1 = psi + dt*((Vo./xLf).*exp(-psi) - abs(Vf)./xLf)

    elseif IDstate == 2
        VdtL = abs(Vf)*dt/xLf
        if VdtL < 1e-6
            psi1 = log( exp(psi-VdtL) + Vo*dt/xLf -
                        0.5*Vo*abs(Vf)*dt*dt/(xLf^2))
        else
            psi1 = log(exp(psi-VdtL) + (Vo/abs(Vf))*(1-exp(-VdtL)))
        end

    elseif IDstate == 3
        psi1 = exp(-abs(Vf)*dt/xLf) * log(abs(Vf)/Vo) + 
        exp(-abs(Vf)*dt/xLf)*psi + log(Vo/abs(Vf))

        if ~any(imag(psi1)) == 0
            return
        end
    end

    return psi1

end
