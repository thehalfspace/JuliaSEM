# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end


# Compute rate-state friciton with depth
function fricDepth(cca, ccb, distN, FltX)
    
    # Friction with depth
    #a_b = cca - ccb

    amax = 0.025; h1 = -15e3/distN
    ao = 0.010; h2 = -18e3/distN

    fP1 = [ao, h1]
    fP2 = [amax, h2]

    fric_depth1 = find(abs.(FltX) .<= 15e3)
    fric_depth2 = find(abs(fP1[2]) .< abs.(FltX) .<= abs(fP2[2]))
    fric_depth3 = find(abs.(FltX) .>= 18e3)
    

    cca[fric_depth1] =  ao   
    cca[fric_depth2] = Int1D(fP1, fP2, FltX[fric_depth2])
    cca[fric_depth3] = amax

    return cca, ccb
end

# Shear stress
function tauDepth(distN, FltX, Vinit, Vo, fo, Seff, cca, ccb)
   
    tmp1 = exp((fo + ccb*log(Vo/Vinit))/cca) 
    
    tauo = Seff*cca*asinh((Vinit/(2*Vo))*tmp1) + Vinit

    return tauo
end

