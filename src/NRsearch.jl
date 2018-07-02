####################################
#   NEWTON RHAPSON SEARCH METHOD
####################################

function NRsearch(fo, Vo, cca, ccb, Seff, tau, tauo, psi, FltZ, FltVfree)

    Vw = 1e10
    fact = 1 + (Vo/Vw)*exp(-psi)

    # NR search point by point for tau if Vf < Vlimit
    eps = 0.001*cca*Seff
    k = 0
    delta = Inf

    while abs(delta) > eps
        fa = fact*tau/(Seff*cca)
        help = -(fo + ccb*psi)/cca

        help1 = exp(help + fa)
        help2 = exp(help - fa)

        Vf = Vo*(help1 - help2)

        Vfprime = fact*(Vo/(cca*Seff))*(help1 + help2)

        delta = (FltZ*FltVfree - FltZ*Vf + tauo - tau)/(1 + FltZ*Vfprime)

        tau = tau + delta
        k = k + 1

        if abs(delta) > 1e10 || k == 1000
            println(k)
            error("NR search fails to converge")
        end
    end

    fa = fact*tau/(Seff*cca)
    
    help = -(fo + ccb*psi)/cca

    help1 = exp(help + fa)
    help2 = exp(help - fa)

    Vf = Vo*(help1 - help2)

    return Vf, tau
end