#######################################################################
#	PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION
#######################################################################

using Parameters

@with_kw struct parameters
    
    # Domain size
    Nsize::Int = 2

    LX::Int = 24e3*Nsize  # depth dimension of rectangular domain
    LY::Int = 15e3*Nsize # off fault dimenstion of rectangular domain

    NelX::Int = 15*Nsize # no. of elements in x
    NelY::Int = 10*Nsize # no. of elements in y

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
    
    Total_time::Int = 100*yr2sec     # Set the total time for simulation here

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

    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    mu::Matrix{Float64} = zeros(NGLL, NGLL)


    # Low velocity layer dimensions
    ThickX::Float64 = LX - LX#8e3
    ThickY::Float64 = 0.75e3

    #.......................
    # EARTHQUAKE PARAMETERS
    #.......................

    Vpl::Float64 = 35e-3/yr2sec	#	Plate loading

    fo::Array{Float64} 	= repeat([0.6], FltNglob)		#	Reference friction coefficient
    Vo::Array{Float64} 	= repeat([1e-6], FltNglob)		#	Reference velocity 'Vo'
    xLf::Array{Float64} = repeat([0.008], FltNglob)#	Dc (Lc) = 8 mm

    Vthres::Float64 = 0.01
    Vevne::Float64 = Vthres

end
