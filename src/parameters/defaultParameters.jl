#######################################################################
#	PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION
#######################################################################


struct parameters

    # Domain size
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

    # Jacobian for global -> local coordinate conversion
    dx_dxi::Float64
    dy_deta::Float64
    jac::Float64
    coefint1::Float64
    coefint2::Float64

    #..................
    # TIME PARAMETERS
    #..................
    yr2sec::Int
    Total_time::Int
    CFL::Float64
    IDstate::Int

    # Some other time variables used in the loop
    dtincf::Float64
    gamma_::Float64
    tevneinc::Int
    dtmax::Int

    #...................
    # MEDIUM PROPERTIES
    #...................

    rho1::Float64
    vs1::Float64

    rho2::Float64
    vs2::Float64

    ETA::Float64

    # Low velocity layer dimensions
    ThickX::Float64
    ThickY::Float64

    #.......................
    # EARTHQUAKE PARAMETERS
    #.......................

    Vpl::Float64

    fo::Array{Float64}
    Vo::Array{Float64}
    xLf::Array{Float64}

    Vthres::Float64
    Vevne::Float64

end



function setParameters(FZdepth)
    
    # Domain size
    Nsize::Int = 2

    LX::Int = 24e3*Nsize  # depth dimension of rectangular domain
    LY::Int = 15e3*Nsize # off fault dimenstion of rectangular domain

    NelX::Int = 60*Nsize # no. of elements in x
    NelY::Int = 40*Nsize # no. of elements in y

    dxe::Float64 = LX/NelX #	Size of one element along X
    dye::Float64 = LY/NelY #	Size of one element along Y
    Nel::Int = NelX*NelY # Total no. of elements
    
    P::Int = 4		#	Lagrange polynomial degree
    NGLL::Int = P + 1 #	No. of Gauss-Legendre-Lobatto nodes
    FltNglob::Int = NelX*(NGLL - 1) + 1

    # Jacobian for global -> local coordinate conversion
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    coefint1::Float64 = jac/dx_dxi^2
    coefint2::Float64 = jac/dy_deta^2

    #..................
    # TIME PARAMETERS
    #..................

    yr2sec::Int = 365*24*60*60
    
    Total_time::Int = 200*yr2sec     # Set the total time for simulation here

    CFL::Float64 = 0.6	#	Courant stability number
     
    IDstate::Int = 2    #   State variable equation type

    # Some other time variables used in the loop
    dtincf::Float64 = 1.2
    gamma_::Float64 = pi/4
    tevneinc::Int = 1
    dtmax::Int = 100 * 24 * 60*60		# 100 days


    #...................
    # MEDIUM PROPERTIES
    #...................

    rho1::Float64 = 2670
    vs1::Float64 = 3464

    rho2::Float64 = 2500
    vs2::Float64 = 0.6*vs1

    ETA = 0



    # Low velocity layer dimensions
    ThickX::Float64 = LX - ceil(FZdepth/dxe)*dxe # ~FZdepth m deep
    ThickY::Float64 = ceil(0.5e3/dye)*dye   # ~ 0.75*2 km wide

    #.......................
    # EARTHQUAKE PARAMETERS
    #.......................

    Vpl::Float64 = 35e-3/yr2sec	#	Plate loading

    fo::Array{Float64} 	= repeat([0.6], FltNglob)		#	Reference friction coefficient
    Vo::Array{Float64} 	= repeat([1e-6], FltNglob)		#	Reference velocity 'Vo'
    xLf::Array{Float64} = repeat([0.008], FltNglob)#	Dc (Lc) = 8 mm

    Vthres::Float64 = 0.01
    Vevne::Float64 = Vthres


    return parameters(Nsize, LX, LY, NelX, NelY, dxe, dye, Nel, P, NGLL, 
                      FltNglob, dx_dxi, dy_deta, jac, coefint1, coefint2, 
                      yr2sec, Total_time, CFL, IDstate, dtincf, gamma_,
                      tevneinc, dtmax, rho1, vs1, rho2, vs2, ETA, ThickX,
                      ThickY, Vpl, fo, Vo, xLf, Vthres, Vevne)

end
