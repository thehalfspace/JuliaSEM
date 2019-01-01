# Calculate XiLf used in computing the timestep
function XiLfFunc(P::params, muMax, cca, ccb, Seff)

    hcell = P.LX/(P.FltNglob-1)
    Ximax = 0.5
    Xithf = 1

    Xith:: Array{Float64} = zeros(P.FltNglob)
    XiLf::Array{Float64} = zeros(P.FltNglob)

    @inbounds for j = 1:P.FltNglob

        # Compute time restricting params
        expr1 = -(cca[j] - ccb[j])/cca[j]
        expr2 = P.gamma_*muMax/hcell*P.xLf[j]/(cca[j]*Seff[j])
        ro = expr2 - expr1

        if (0.25*ro*ro - expr2) >= 0
            Xith[j] = 1/ro
        else
            Xith[j] = 1 - expr1/expr2
        end

        # For each node, compute slip that node cannot exceed in one timestep
        if Xithf*Xith[j] > Ximax
            XiLf[j] = Ximax*P.xLf[j]
        else
            XiLf[j] = Xithf*Xith[j]*P.xLf[j]
        end
    end
       

    return XiLf
end


# K diagonal vector computation
function KdiagFunc(P::params, iglob, W, H, Ht, FltNI)

	nglob = P.FltNglob*(P.NelY*(P.NGLL-1) + 1)

    # Compute the diagonal of K
    Kdiag::Array{Float64} = zeros(nglob)
    Klocdiag::Array{Float64,2} = zeros(P.NGLL, P.NGLL)
    @inbounds for et = 1:P.Nel
        ig = iglob[:,:,et]
        wloc = W[:,:,et]
        Klocdiag[:,:] .= 0

        for k =  1:P.NGLL
            for j = 1:P.NGLL
                Klocdiag[k,j] +=  sum( P.coefint1*H[k,:].*(wloc[:,j].*Ht[:,k])
                                + P.coefint2*(wloc[k,:].*H[j,:]).*Ht[:,j] )
            end
        end

        Kdiag[ig] .+= Klocdiag[:,:]
    end

    return Kdiag[FltNI]

end

# IDstate functions
function IDS(xLf, Vo, psi, dt, Vf, cnd, IDstate = 2)
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
function IDS2(xLf, Vo, psi, psi1, dt, Vf, Vf1, IDstate = 2)
            
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

# Slip rates on fault for quasi-static regime
function slrFunc(P, NFBC, FltNglob, psi, psi1, Vf, Vf1, IDstate, tau1, tauo, Seff, cca, ccb, dt)


    tauAB = SharedArray{Float64}(FltNglob)

   @sync @distributed for j = NFBC:FltNglob 

        psi1[j] = IDS(P.xLf[j], P.Vo[j], psi[j], dt, Vf[j], 1e-6, IDstate)

        tauAB[j] = tau1[j] + tauo[j]
        fa = tauAB[j]/(Seff[j]*cca[j])
        help = -(P.fo[j] + ccb[j]*psi1[j])/cca[j]
        help1 = exp(help + fa)
        help2 = exp(help - fa)
        Vf1[j] = P.Vo[j]*(help1 - help2) 
    end

    return psi1, Vf1
end
