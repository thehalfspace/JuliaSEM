
function main(P::parameters, S::input_variables)

    Ht = S.H'
    wgll2 = S.wgll*S.wgll';
    
    # Time solver variables
    dt = S.dt0
    dtmin = dt
    half_dt = 0.5*dtmin
    half_dt_sq = 0.5*dtmin^2

    # dt modified slightly for damping
    if P.ETA != 0
	    dt = dt/sqrt(1 + 2*P.ETA)
    end

    # Initialize kinematic field: global arrays
    global d = zeros(S.nglob)
    global v = zeros(S.nglob)
    v .= 0.5e-3
    global a = zeros(S.nglob)

    #.....................................
    # Stresses and time related variables
    #.....................................
    tau::Array{Float64} = zeros(P.FltNglob)
    FaultC::Array{Float64} = zeros(P.FltNglob)
    Vf1::Array{Float64}  = zeros(P.FltNglob)
    Vf2::Array{Float64} = zeros(P.FltNglob)
    Vf0::Array{Float64} = zeros(length(S.iFlt))
    FltVfree::Array{Float64} = zeros(length(S.iFlt))
    psi::Array{Float64} = zeros(P.FltNglob)
    psi0::Array{Float64} = zeros(P.FltNglob)
    psi1::Array{Float64} = zeros(P.FltNglob)
    psi2::Array{Float64} = zeros(P.FltNglob)
    tau1::Array{Float64} = zeros(P.FltNglob)
    tau2::Array{Float64} = zeros(P.FltNglob)
    tau3::Array{Float64} = zeros(P.FltNglob)
    tauNR::Array{Float64} = zeros(P.FltNglob)

    # Initial state variable
    psi = S.tauo./(S.Seff.*S.ccb) - P.fo./S.ccb - (S.cca./S.ccb).*log.(2*v[S.iFlt]./P.Vo)
    psi0 .= psi[:]

    isolver::Int = 1
    
    # Some more initializations
    r::Array{Float64} = zeros(S.nglob)
    beta_::Array{Float64} = zeros(S.nglob)
    alpha_::Array{Float64} = zeros(S.nglob)

    F::Array{Float64} = zeros(S.nglob)
    dPre::Array{Float64} = zeros(S.nglob)
    vPre::Array{Float64} = zeros(S.nglob)
    dd::Array{Float64} = zeros(S.nglob)
    dacum::Array{Float64} = zeros(S.nglob)
    dnew::Array{Float64} = zeros(length(S.FltNI))

    # Iterators
    idelevne = 3
    tevneb = 0
    tevne = 0

    v = v[:] .- 0.5*P.Vpl
    Vf::Array{Float64} = 2*v[S.iFlt]
    iFBC::Array{Int64} = findall(abs.(S.FltX) .> 24e3)
    NFBC::Int64 = length(iFBC)
    Vf[iFBC] .= 0


    v[S.FltIglobBC] .= 0

    # Preallocate variables with unknown size
    time_ = zeros(1000000)

    Stress::Array{Float64} = zeros(P.FltNglob, 1000000)
    SlipVel::Array{Float64} = zeros(P.FltNglob, 1000000)
    Slip::Array{Float64} = zeros(P.FltNglob, 1000000)


    #....................
    # Start of time loop
    #....................
    it = 0
    t = 0

    it = it + 1
    t = t + dt

    time_[it] = t 

    vPre .= v
    dPre .= d

    Vf0 .= 2*v[S.iFlt] .+ P.Vpl
    Vf  .= Vf0
        
    # Compute the forcing term
    F .= 0
    F[S.iFlt] .= dPre[S.iFlt] .+ v[S.iFlt]*dt

    # Assign previous displacement field as initial guess
    dnew .= d[S.FltNI]

    # Solve d = K^-1F by PCG
    #println("\nPCG:")
    dnew = PCG(P, S.diagKnew, dnew, F, S.iFlt, S.FltNI,
                  S.H, Ht, S.iglob, S.nglob, S.W)

    return dnew
end
