#########################################
#										#
#	SPECTRAL ELEMENT METHOD FOR  		#
#	EARTHQUAKE CYCLE SIMULATION			#
#										#
#	Prithvi Thakur						#
#	Adapted from Kaneko et al. (2011)	#
#	and J.P. Ampuero matlab exercises	#
#########################################

#.................................
# Include external function files
#.................................
include("Parameters.jl")	#	Set Parameters
include("GetGLL.jl")		#	Polynomial interpolation
include("Meshbox.jl")		# 	Build 2D mesh
include("BoundaryMatrix.jl") #	Boundary matrices
include("FindNearestNode.jl")#	Nearest node	


#...............
# Build 2D Mesh
#...............
iglob, x, y = MeshBox(LX, LY, NelX, NelY, NGLL)
x = x-LX		#	For halfspace
nglob = length(x);


#..................
# Initialization
#..................
# The derivatives of the Lagrange Polynomials were pre-tabulated 
# xgll = location of the GLL nodes inside the reference segment [-1,1]
xgll, wgll, H = GetGLL(NGLL)
Ht = H'
wgll2 = wgll*wgll';

# For internal forces
W = zeros(NGLL, NGLL, Nel)

# Global Mass Matrix
M = zeros(nglob, 1)

# Mass+Damping matrix
MC = zeros(nglob, 1);

# Used for variable time stepping
muMax = 0

# Assemble the Mass and Stiffness matrices
include("Assemble.jl")

# Time solver parameters
dt = CFL*dt
half_dt = 0.5*dt
dtmin = dt

tmas = Total_time
dtmax = 5 * 24 * 60*60/distN * 1000		# 5 days

# dt modified slightly for damping
if ETA != 0
	dt = dt/sqrt(1 + 2*ETA)
end

# Initialize kinematic field: global arrays
d = zeros(nglob, 1)
v = zeros(nglob, 1)
v[:,:] = 0.5e-3
a = zeros(nglob,1)



#......................
# Boundary conditions : absorbing boundaries on 3 sides, fault boundary on one side
#......................

# Left boundary
impedance = rho1*vs1
BcLC, iBcL = BoundaryMatrix(wgll, NelX, NelY, iglob, dy_deta, impedance, 'L')

# Right Boundary = free surface: nothing to do

# Top Boundary
impedance = rho1*vs1
BcTC, iBcT = BoundaryMatrix(wgll, NelX, NelY, iglob, dx_dxi, impedance, 'T')

# Mass matrix at boundaries
Mq = M
M[iBcL] = M[iBcL] + half_dt*BcLC
M[iBcT] = M[iBcT] + half_dt*BcTC


# Dynamic fault at bottom boundary
FltB, iFlt = BoundaryMatrix(wgll, NelX, NelY, iglob, dx_dxi, 1, 'B') # impedance = 1
FltNglob = NelX*(NGLL - 1) + 1

FltZ = M[iFlt]./FltB /half_dt * 0.5
FltX = x[iFlt]



# Linear interpolation function
function Int1D(P1, P2, val)	
	Line = P1[1] .+ ( (P2[1] - P1[1])/((P2[2] - P1[2])).*(val .- P1[2]) )	
	return Line
end


#.....................
# Initial Conditions 
#.....................
tauoBack = 22.5e6
tauo = repmat([tauoBack], FaultNglob, 1)


# Friction with depth
a_b = cca - ccb
fP1 = [0 -1.2e3/distN]
fP2 = [-0.0041 -2e3/distN]
fP3 = [-0.0041 -12e3/distN]
fP4 = [0.015 -17e3/distN]
fP5 = [0.024 -24e3/distN]

fric_depth1 = find(abs.(FltX) .<= abs(fP2[2]))
fric_depth2 = find(abs(fP2[2]) .< abs.(FltX) .<= abs(fP3[2]))
fric_depth3 = find(abs(fP3[2]) .< abs.(FltX) .<= abs(fP4[2]))
fric_depth4 = find(abs(fP4[2]) .< abs.(FltX) .<= abs(fP5[2]))
fric_depth5 = find(abs.(FltX) .> abs(fP5[2]))

a_b[fric_depth1] = Int1D(fP1, fP2, FltX[fric_depth1])
a_b[fric_depth2] = Int1D(fP2, fP3, FltX[fric_depth2])
a_b[fric_depth3] = Int1D(fP3, fP4, FltX[fric_depth3])
a_b[fric_depth4] = Int1D(fP4, fP5, FltX[fric_depth4])
a_b[fric_depth5] = 0.0047

cca = ccb + a_b

# Effective normal stress
sP1 = [3e6 0]
sP2 = [50e6 -2e3/distN]
Seff_depth = find(abs.(FltX) .<= abs(sP2[2]))
Seff[Seff_depth] = Int1D(sP1, sP2, FltX[Seff_depth])


# Initial shear stress
tP1 = [3e6 0]
tP2 = [30e6 -2e3/distN]
tP3 = [30e6 -12e3/distN]
tP4 = [22.5e6 -17e3/distN]
tP5 = [22.5e6 -24e3distN]

tau_depth1 = find(abs.(FltX) .<= abs(tP2[2]))
tau_depth2 = find(abs(tP2[2]) .< abs.(FltX) .<= abs(tP3[2]))
tau_depth3 = find(abs(tP3[2]) .< abs.(FltX) .<= abs(tP4[2]))
tau_depth4 = find(abs(tP4[2]) .< abs.(FltX) .<= abs(tP5[2]))

tauo[tau_depth1] = Int1D(tP1, tP2, FltX[tau_depth1])
tauo[tau_depth2] = Int1D(tP2, tP3, FltX[tau_depth2])
tauo[tau_depth3] = Int1D(tP3, tP4, FltX[tau_depth3])
tauo[tau_depth4] = Int1D(tP4, tP5, FltX[tau_depth4])



#.....................................
# Stresses and time related variables
#.....................................
tau = zeros(FaultNglob, 1)

psi = tauo./(Seff.*ccb) - fo./ccb - (cca./ccb).*log.(2*v[iFlt]./Vo)
psi0 = psi

# Kelvin-Voigt Viscosity
if ETA !=0
	Nel_ETA = NelX
	x1 = 0.5*(1 + xgll')
	eta_taper = exp.(-pi*x1.^2)
	eta = ETA*dt * repmat([eta_taper], NGLL, 1)
end


# Output seismograms
OUTxseis, OUTyseis, OUTiglob, OUTdseis = FindNearestNode(OUTxseis, OUTyseis, x, y)
kkseis = 1

# Initialize data for output snapshots
#OUTindx = 

# Time related variables
t = 0
it = 0
dtincf = 1.2
dtpre = dt
gamma_ = pi/4

# Avg. node spacing
hcell = LX/(FaultNglob - 1)
Ximax = 0.5
Xithf = 1
trec = 0
Vthres =  0.01
slipstart = 0
ievb = 0
ieva = 0
ntvsx = 0
tvsx = 0.2*yr2sec #	Output slip every 0.2 years
tvsxinc = tvsx
nevne = 0
tevneinc = 0.5
Vevne = Vthres

# Compute XiLf for each fault node
Xith = zeros(FaultNglob,1)
XiLf = zeros(FaultNglob,1)

for j=1:FaultNglob
	
	# Compute time restricting parameters
	expr1 = -(cca[j] - ccb[j])/cca[j]
	expr2 = gamma_*muMax/hcell*xLf[j]/(cca[j]*Seff[j])
	ro = expr2 - expr1

	if (0.25*ro*ro-expr2) >= 0
		Xith[j] = 1/ro
	else
		Xith[j] = 1 - expr1/expr2
	end

	# For each node, compute slip that node cannot exceed in one timestep
	if Xithf*Xith[j] > Ximax
		XiLf[j] = Ximax*xLf[j]
	else
		XiLf[j] = Xithf*Xith[j]*xLf[j]
	end
end

# Time related variables added
OUTt = 0.5
q = 1
OUTtGo = 0
OUTtCount = 1

# Output variables at several locations on fault
OutLoc1 = 1e3/distN
OutLoc2 = 2e3/distN
OutLoc3 = 6e3/distN
OutLoc4 = 8e3/distN

xLoc1, yLoc1, iglobLoc1 = FindNearestNode(OutLoc1, 0, x, y)
FltLoc1 = round((LX - OutLoc1)*FaultNglob/LX)

xLoc1, yLoc1, iglobLoc2 = FindNearestNode(OutLoc2, 0, x, y)
FltLoc2 = round((LX - OutLoc2)*FaultNglob/LX)

xLoc1, yLoc1, iglobLoc3 = FindNearestNode(OutLoc3, 0, x, y)
FltLoc3 = round((LX - OutLoc3)*FaultNglob/LX)

xLoc1, yLoc1, iglobLoc4 = FindNearestNode(OutLoc4, 0, x, y)
FltLoc4 = round((LX - OutLoc4)*FaultNglob/LX)

# Display important parameters 
println("Total number of nodes on fault = ", FaultNglob)
println("Average node spacing = ", LX/distN/(FaultNglob-1))
@printf("dt = %1.17f", dt)


# Find nodes that do not belong to fault
FltNI = deleteat!(collect(1:nglob), iFlt) 

# Some more initializations
r = zeros(nglob,1)
beta_ = zeros(nglob,1)
alpha_ = zeros(nglob,1)
p = zeros(nglob, 1)

F = zeros(nglob,1)
dPre = zeros(nglob,1)
vPre = zeros(nglob,1)
dd = zeros(nglob,1)
dacum = zeros(nglob,1)

dnew = zeros(length(FltNI),1)

println("\nSetup Complete")
