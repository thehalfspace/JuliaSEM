#########################
# SOME OTHER FUNCTIONS
#########################

function IDS(psi, dt, Vo, xLf, Vf, cnd, IDstate = 2)
    #= compute slip-rates on fault based on different
       formulations =#

    if IDstate == 1
        psi1 = psi + dt*((Vo./xLf).*exp(-psi) - abs(Vf)./xLf)

    elseif IDstate == 2
        VdtL = abs(Vf)*dt/xLf
        if VdtL < cnd
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

# On fault slip rates
function IDS2(psi, psi1, dt, Vo, xLf, Vf, Vf1, IDstate = 2)
            
    if IDstate == 1
        psi2 = psi + 0.5*dt*( (Vo/xLf)*exp(-psi) - abs(Vf)/xLf 
                                + (Vo/xLf)*exp(-psi1) - abs(Vf1)/xLf )

    elseif IDstate == 2
        VdtL = 0.5*abs(Vf1 + Vf)*dt/xLf

        if VdtL < 1e-6
            psi2 = log( exp(psi-VdtL) + Vo*dt/xLf -
                            0.5*Vo*0.5*abs(Vf1 + Vf)*dt*dt/(xLf^2))
        else
            psi2 = log(exp(psi-VdtL) + 
                            (Vo/(0.5*abs(Vf + Vf1)))*(1-exp(-VdtL)))
        end

    elseif IDstate == 3
        psi2 = exp(-0.5*abs(Vf + Vf1)*dt/xLf) * log(0.5*abs(Vf + Vf1)/Vo) + 
                        exp(-0.5*abs(Vf + Vf1)*dt/xLf)*psi 
                        + log(Vo/(-0.5*abs(Vf + Vf1)) )
    end

    return psi2
end

