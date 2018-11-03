#################################
# READ OUTPUT FROM SIMULATION
#################################

mutable struct results
    Stress::Array{Float64,2}
    SlipVel::Array{Float64,2}
    Slip::Array{Float64,2}
    time_::Array{Float64}
end

struct parameters
    
    Nsize::Int 
    LX::Int
    LY::Int
    NelX::Int
    NelY::Int

    dxe::Float64
    dye::Float64
    Nel::Int
    
    P::Int    
    NGLL::Int
    FltNglob::Int

    dx_dxi::Float64
    dy_deta::Float64
    jac::Float64
    coefint1::Float64
    coefint2::Float64

    yr2sec::Int
    Total_time::Int
    CFL::Float64
    IDstate::Int

    dtincf::Float64
    gamma_::Float64
    tevneinc::Int
    dtmax::Int

    rho1::Float64
    vs1::Float64

    rho2::Float64
    vs2::Float64

    ETA::Float64

    ThickX::Float64
    ThickY::Float64

    Vpl::Float64

    fo::Array{Float64}
    Vo::Array{Float64}
    xLf::Array{Float64}

    Vthres::Float64
    Vevne::Float64

end

struct input_variables
    
    iglob::Array{Int,3}
    nglob::Int
    x::Array{Float64}
    y::Array{Float64}

    xgll::Array{Float64}
    wgll::Array{Float64}
    H::Array{Float64,2}
    Ht::Array{Float64,2}

    W::Array{Float64,3}
    M::Array{Float64}
    MC::Array{Float64}

    BcLC::Array{Float64}
    iBcL::Array{Int}
    BcTC::Array{Float64}
    iBcT::Array{Int}

    FltB::Array{Float64}
    iFlt::Array{Int64}
    FltZ::Array{Float64}
    FltX::Array{Float64}

    cca::Array{Float64}
    ccb::Array{Float64}
    Seff::Array{Float64}
    tauo::Array{Float64}

    Nel_ETA::Float64
    XiLf::Array{Float64}
    diagKnew::Array{Float64}

    FltIglobBC::Array{Int}
    FltNI::Array{Int}

    dt0::Float64

end

using Serialization
open("output/test12.out") do f
    global O, sim_time, P, S
    O = deserialize(f)
    sim_time = deserialize(f)
    P = deserialize(f)
    S = deserialize(f)
end
