########################################################
#
#	FUNCTION TO BUILD THE 2D BOUNDARIES FOR SEM
#
########################################################

function BoundaryMatrix(s::space_parameters, m::medium_properties, wgll, iglob, side)

	# INPUT: 
	#		wgll = GLL weights (see GetGLL)
	#		NelX, NelY = no. of elements
	#		iglob = global to local index
	#		jac1D = line Jacobian
	#		side = 'L', 'R', 'T', 'B'


	if side == 'L'
		eB = collect(0:s.NelY-1)*s.NelX .+ 1
		igll = 1
		jgll = collect(1:s.NGLL)
        jac1D = s.dy_deta
        impedance = m.rho1*m.vs1

	elseif side == 'R'
		eB = collect(0:s.NelY-1)*s.NelX .+ s.NelX
		igll = NGLL
		jgll = collect(1:s.NGLL)
        jac1D = s.dy_deta
        impedance = m.rho1*m.vs1

	elseif side == 'T'
		eB = (s.NelY-1)*s.NelX .+ collect(1:s.NelX)
		igll = collect(1:s.NGLL)
		jgll = s.NGLL
        jac1D = s.dx_dxi
        impedance = m.rho1*m.vs1

	else 
		eB = collect(1:s.NelX)
		igll = collect(1:s.NGLL)
		jgll = 1
        jac1D = s.dx_dxi
        impedance = 1
	end

	NelB = length(eB)
	ng = NelB*(s.NGLL-1) .+ 1
	iB = zeros(Int32, ng)
	B = zeros(ng)
	jB = zeros(s.NGLL, NelB)

	for e=1:NelB
		ip = (s.NGLL-1)*(e-1) .+ collect(1:s.NGLL)
		iB[ip] = iglob[igll, jgll, eB[e]]
		jB[:,e] = ip
		B[ip] .+= jac1D*wgll*impedance
	end

	return B, iB
end
