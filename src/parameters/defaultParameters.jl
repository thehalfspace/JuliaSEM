#######################################################################
#	PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION
#######################################################################

using Parameters

@with_kw struct space_parameters
    
    # Domain size
    Nsize::Int8 = 2

    LX::Int64 = 24e3*Nsize  # depth dimension of rectangular domain
    LY::Int64 = 15e3*Nsize # off fault dimenstion of rectangular domain

    NelX::Int8 = 15*Nsize # no. of elements in x
    NelY::Int8 = 10*Nsize # no. of elements in y

    dxe::Float64 = LX/NelX #	Size of one element along X
    dye::Float64 = LY/NelY #	Size of one element along Y
    Nel::Int64 = NelX*NelY # Total no. of elements
    
    P::Int = 4		#	Lagrange polynomial degree
    NGLL::Int64 = P + 1 #	No. of Gauss-Legendre-Lobatto nodes
    FltNglob::Int64 = NelX*(NGLL - 1) + 1

    # Jacobian for global -> local coordinate conversion
    dx_dxi::Float64 = 0.5*dxe
    dy_deta::Float64 = 0.5*dye
    jac::Float64 = dx_dxi*dy_deta
    coefint1::Float64 = jac/dx_dxi^2
    coefint2::Float64 = jac/dy_deta^2
end

@with_kw struct time_parameters
    
    yr2sec::Int64 = 365*24*60*60
    Total_time::Int128 = 60*yr2sec     # Set the total time for simulation here

    CFL::Float64 = 0.6	#	Courant stability number
     
    IDstate::Int = 2    #   State variable equation type

    # Some other time variables used in the loop
    dtincf::Float64 = 1.2
    gamma_::Float64 = pi/4
    tevneinc::Int64 = 1
    dtmax::Int64 = 100 * 24 * 60*60		# 100 days

end



@with_kw struct medium_properties

    NGLL = space_parameters().NGLL
    LX = space_parameters().LX

    rho1::Float64 = 2670
    vs1::Float64 = 3464

    rho2::Float64 = 2500
    vs2::Float64 = 0.6*vs1

    ETA = 0

    rho::Matrix{Float64} = zeros(NGLL, NGLL)
    mu::Matrix{Float64} = zeros(NGLL, NGLL)


    # Low velocity layer dimensions
    ThickX::Float64 = LX - 8e3
    ThickY::Float64 = 0.75e3
end


@with_kw struct earthquake_parameters

    yr2sec = time_parameters().yr2sec
    FltNglob = space_parameters().FltNglob

    Vpl::Float64 = 35e-3/yr2sec	#	Plate loading

    fo::Array{Float64} 	= repeat([0.6], FltNglob)		#	Reference friction coefficient
    Vo::Array{Float64} 	= repeat([1e-6], FltNglob)		#	Reference velocity 'Vo'
    xLf::Array{Float64} = repeat([0.008], FltNglob)#	Dc (Lc) = 8 mm

    Vthres::Float64 = 0.01
    Vevne::Float64 = Vthres

end
