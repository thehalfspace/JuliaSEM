########################################################
#
#	FUNCTION TO BUILD THE 2D BOUNDARIES FOR SEM
#
########################################################

function BoundaryMatrix(wgll, NelX, NelY, iglob, jac1D, impedance, side)

	# INPUT: 
	#		wgll = GLL weights (see GetGLL)
	#		NelX, NelY = no. of elements
	#		iglob = global to local index
	#		jac1D = line Jacobian
	#		side = 'L', 'R', 'T', 'B'

	NGLL = length(wgll)

	if side == 'L'
		eB = collect(0:NelY-1)*NelX + 1
		igll = 1
		jgll = collect(1:NGLL)

	elseif side == 'R'
		eB = collect(0:NelY-1)*NelX + NelX
		igll = NGLL
		jgll = collect(1:NGLL)

	elseif side == 'T'
		eB = (NelY-1)*NelX + collect(1:NelX)
		igll = collect(1:NGLL)
		jgll = NGLL

	else 
		eB = collect(1:NelX)
		igll = collect(1:NGLL)
		jgll = 1
	end

	NelB = length(eB)
	ng = NelB*(NGLL-1) + 1
	iB = zeros(Int32, ng)
	B = zeros(ng)
	jB = zeros(NGLL, NelB)

	for e=1:NelB
		ip = (NGLL-1)*(e-1) + collect(1:NGLL)
		iB[ip] = iglob[igll, jgll, eB[e]]
		jB[:,e] = ip
		B[ip] = B[ip] + jac1D*wgll*impedance
	end

	return B, iB
end
