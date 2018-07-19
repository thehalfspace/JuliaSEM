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


function MeshBox(LX, LY, NELX, NELY, NGLL, dxe, dye) 

	XGLL = GetGLL(NGLL)[1]

	iglob = zeros(Int32, NGLL, NGLL, NELX*NELY)
	nglob = (NELX*(NGLL-1) + 1) * (NELY*(NGLL-1) + 1)

	x = zeros(nglob, 1)
	y = zeros(nglob, 1)

	e = 0
	last_iglob = 0

	ig = reshape(collect(1:NGLL*NGLL), NGLL, NGLL)
	igL = reshape(collect(1:NGLL*(NGLL-1)), NGLL-1, NGLL) # Left edge
	igB = reshape(collect(1:NGLL*(NGLL-1)), NGLL, NGLL-1) # Bottom edge
	igLB = reshape(collect(1:(NGLL-1)*(NGLL-1)), NGLL-1, NGLL-1) # rest of the elements

	x2gll = repmat(0.5*(1+XGLL), 1, NGLL)
	#ygll = dye*xgll'
	#xgll = dxe*xgll


	for ey = 1:NELY         # number of x elements
		for ex = 1:NELX     # number of y elements

			e = e + 1
            
            ygll = dye[ey]*x2gll'
            xgll = dxe[ex]*x2gll

			# Redundant nodes at element edges
            # NGLL = number of GLL nodes per element
            
			if e == 1
				ig = reshape(collect(1:NGLL*NGLL), NGLL, NGLL)
			else
				if ey ==1 # Bottom Row
					ig[1,:] = iglob[NGLL, :, e-1]   # Left edge
					ig[2:end, :] = last_iglob + igL # The rest

				elseif ex == 1	# Left Column
					ig[:,1] = iglob[:,NGLL,e-NELX]	# Bottom edge
					ig[:,2:end] = last_iglob + igB 	# The rest

				else 			# Other Elements
					ig[1,:] = iglob[NGLL, :, e-1]	# Left edge
					ig[:,1] = iglob[:, NGLL, e-NELX]# Bottom edge
					ig[2:end, 2:end] = last_iglob + igLB
				end
			end

			iglob[:,:,e] = ig
			last_iglob = ig[NGLL, NGLL]

			# Global coordinates of computational nodes
            #x[ig] = dxe[ex]*(ex-1) + xgll
            #y[ig] = dye[ey]*(ey-1) + ygll
            
            x[ig] = LX*sin(pi/2*(ex-1)/NelX) + xgll
            y[ig] = LY*(1 - sin(pi/2*(ey-1)/NelY)) + ygll

		end
	end


	return iglob, x, y

end

