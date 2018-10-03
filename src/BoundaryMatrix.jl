########################################################
#
#	FUNCTION TO BUILD THE 2D BOUNDARIES FOR SEM
#
########################################################

function BoundaryMatrix(P::parameters, wgll, iglob, side)

	# INPUT: 
	#		wgll = GLL weights (see GetGLL)
	#		NelX, NelY = no. of elements
	#		iglob = global to local index
	#		jac1D = line Jacobian
	#		side = 'L', 'R', 'T', 'B'


	if side == 'L'
		eB = collect(0:P.NelY-1)*P.NelX .+ 1
		igll = 1
		jgll = collect(1:P.NGLL)
        jac1D = P.dy_deta
        impedance = P.rho1*P.vs1

	elseif side == 'R'
		eB = collect(0:P.NelY-1)*P.NelX .+ P.NelX
		igll = P.NGLL
		jgll = collect(1:P.NGLL)
        jac1D = P.dy_deta
        impedance = P.rho1*P.vs1

	elseif side == 'T'
		eB = (P.NelY-1)*P.NelX .+ collect(1:P.NelX)
		igll = collect(1:P.NGLL)
		jgll = P.NGLL
        jac1D = P.dx_dxi
        impedance = P.rho1*P.vs1

	else 
		eB = collect(1:P.NelX)
		igll = collect(1:P.NGLL)
		jgll = 1
        jac1D = P.dx_dxi
        impedance = 1
	end

	NelB = length(eB)
	ng = NelB*(P.NGLL-1) .+ 1
	iB = zeros(Int32, ng)
	B = zeros(ng)
	jB = zeros(P.NGLL, NelB)

	for e=1:NelB
		ip = (P.NGLL-1)*(e-1) .+ collect(1:P.NGLL)
		iB[ip] = iglob[igll, jgll, eB[e]]
		jB[:,e] = ip
		B[ip] .+= jac1D*wgll*impedance
	end

	return B, iB
end
