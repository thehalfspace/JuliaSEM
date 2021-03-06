# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end


# Compute rate-state friciton with depth
function fricDepth(cca, ccb, distN, FltX)
    
    # Friction with depth
    a_b = cca - ccb
    fP1 = [0 -1.2e3/distN]
    fP2 = [-0.0041 -2e3/distN]
    fP3 = [-0.0041 -12e3/distN]
    fP4 = [0.015 -17e3/distN]
    fP5 = [0.024 -24e3/distN]

    fric_depth1 = find(abs.(FltX) .<= abs(fP2[2]))
    fric_depth2 = find(abs(fP2[2]) .< abs.(FltX) .<= abs(fP3[2]))
    fric_depth3 = find(abs(fP3[2]) .< abs.(FltX) .<= abs(fP4[2]))
    fric_depth4 = find(abs(fP4[2]) .< abs.(FltX) .<= abs(fP5[2]))
    fric_depth5 = find(abs.(FltX) .> abs(fP5[2]))

    a_b[fric_depth1] = Int1D(fP1, fP2, FltX[fric_depth1])
    a_b[fric_depth2] = Int1D(fP2, fP3, FltX[fric_depth2])
    a_b[fric_depth3] = Int1D(fP3, fP4, FltX[fric_depth3])
    a_b[fric_depth4] = Int1D(fP4, fP5, FltX[fric_depth4])
    a_b[fric_depth5] = 0.0047

    cca = ccb + a_b

    return cca, ccb
end



# Effective normal stress
function SeffDepth(distN, FltX)
    sP1 = [3e6 0]
    sP2 = [50e6 -2e3/distN]
    Seff_depth = find(abs.(FltX) .<= abs(sP2[2]))
    Seff[Seff_depth] = Int1D(sP1, sP2, FltX[Seff_depth])

    return Seff
end


# Shear stress
function tauDepth(distN, FltX)

    tP1 = [3e6 0]
    tP2 = [30e6 -2e3/distN]
    tP3 = [30e6 -12e3/distN]
    tP4 = [22.5e6 -17e3/distN]
    tP5 = [22.5e6 -24e3distN]

    tau_depth1 = find(abs.(FltX) .<= abs(tP2[2]))
    tau_depth2 = find(abs(tP2[2]) .< abs.(FltX) .<= abs(tP3[2]))
    tau_depth3 = find(abs(tP3[2]) .< abs.(FltX) .<= abs(tP4[2]))
    tau_depth4 = find(abs(tP4[2]) .< abs.(FltX) .<= abs(tP5[2]))

    tauo[tau_depth1] = Int1D(tP1, tP2, FltX[tau_depth1])
    tauo[tau_depth2] = Int1D(tP2, tP3, FltX[tau_depth2])
    tauo[tau_depth3] = Int1D(tP3, tP4, FltX[tau_depth3])
    tauo[tau_depth4] = Int1D(tP4, tP5, FltX[tau_depth4])

    return tauo
end

