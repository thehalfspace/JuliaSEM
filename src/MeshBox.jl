###############################################################################
#		Spectral Element Mesh for Rectangular Box, with internal 
#		Gauss-Legendre-Lobatto (GLL) dub-grids.
#
#	INPUT:	LX = x-dimension
#			LY = y-dimension
#			NELX = no. of elements in x
#			NELY = no. of elements in y
#			NGLL = no. of GLL nodes


#	OUTPUT:	iglob[ngll, ngll, NELX, NELY] = maps local to global
#											numbering
#			I = iglob[i,j,e] is the global node index of the 
#							 (i,j)th GLL node internal to the
#							 e-th element.

#			Elements are numbered row by row from from bottom-left
#			to top-right. The table iglob is needed to assemble 
#			global data from local data.

#			x[:] = global x coordinates of GLL nodes, starting at 0
#			y[:] = global y coordinates of GLL nodes, starting at 0
###############################################################################


function MeshBox(P::parameters) 

	XGLL = GetGLL(P.NGLL)[1]

	iglob = zeros(Int, P.NGLL, P.NGLL, P.Nel)
	nglob = P.FltNglob*(P.NelY*(P.NGLL-1) + 1)

    x::Array{Float64} = zeros(nglob)
    y::Array{Float64} = zeros(nglob)

	et = 0
	last_iglob = 0

	ig = reshape(collect(1:P.NGLL*P.NGLL), P.NGLL, P.NGLL)
	igL = reshape(collect(1:P.NGLL*(P.NGLL-1)), P.NGLL-1, P.NGLL) # Left edge
	igB = reshape(collect(1:P.NGLL*(P.NGLL-1)), P.NGLL, P.NGLL-1) # Bottom edge
	igLB = reshape(collect(1:(P.NGLL-1)*(P.NGLL-1)), P.NGLL-1, P.NGLL-1) # rest of the elements

	xgll = repeat(0.5*(1 .+ XGLL), 1, P.NGLL)
	ygll = P.dye*xgll'
	xgll = P.dxe*xgll


	for ey = 1:P.NelY         # number of x elements
		for ex = 1:P.NelX     # number of y elements

			et = et + 1

			# Redundant nodes at element edges
            
            # NGLL = number of GLL nodes per element
            
			if et == 1
				ig = reshape(collect(1:P.NGLL*P.NGLL), P.NGLL, P.NGLL)
			else
				if ey ==1 # Bottom Row
					ig[1,:] = iglob[P.NGLL, :, et-1]   # Left edge
					ig[2:end, :] = last_iglob .+ igL # The rest

				elseif ex == 1	# Left Column
					ig[:,1] = iglob[:,P.NGLL,et-P.NelX]	# Bottom edge
					ig[:,2:end] = last_iglob .+ igB 	# The rest

				else 			# Other Elements
					ig[1,:] = iglob[P.NGLL, :, et-1]	# Left edge
					ig[:,1] = iglob[:, P.NGLL, et-P.NelX]# Bottom edge
					ig[2:end, 2:end] = last_iglob .+ igLB
				end
			end

			iglob[:,:,et] = ig
			last_iglob = ig[P.NGLL, P.NGLL]

			# Global coordinates of computational nodes
			x[ig] .= P.dxe*(ex-1) .+ xgll
			y[ig] .= P.dye*(ey-1) .+ ygll

		end
	end

	return iglob, x, y
end

