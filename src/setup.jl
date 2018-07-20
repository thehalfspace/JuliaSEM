#########################################
#										#
#	SPECTRAL ELEMENT METHOD FOR  		#
#	EARTHQUAKE CYCLE SIMULATION			#
#										#
#	Prithvi Thakur						#
#	Adapted from Kaneko et al. (2011)	#
#	and J.P. Ampuero's SEMLAB       	#
#########################################

#.................................
# Include external function files
#.................................
include("parameters/defaultParameters.jl")	#	Set Parameters
include("GetGLL.jl")		#	Polynomial interpolation
include("mesh/VariableMesh.jl")		# 	Build 2D mesh
include("BoundaryMatrix.jl") #	Boundary matrices
include("FindNearestNode.jl")#	Nearest node	
include("IDState.jl") # other functions

#...............
# Build 2D Mesh
#...............
iglob, x, y = MeshBox(LX, LY, NelX, NelY, NGLL, dxe, dye)
x = x-LX		#	For halfspace
x_points = x_points - LX

const nglob = length(x);


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
M = zeros(nglob)

# Mass+Damping matrix
MC = zeros(nglob)

#..................................
# Jacobian with variable mesh size
#..................................
coefint1 = zeros(Nel)
coefint2 = zeros(Nel)

function coefint(coefint1, coefint2, dxe, dye)

    eo = 0
    for ey = 1:length(dye)
        for ex = 1:length(dxe)

            eo = eo + 1

            jac = 0.25*dxe[ex]*ye[ey]

            coefint1[eo] = jac/(0.5*dxe[ex])^2
            coefint2[eo] = jac/(0.5*dye[ey])^2

        end
    end
    
    return coefint1, coefint2
end

coefint1, coefint2 = coefint(coefint1, coefint2, dxe, dye)

# Assemble the Mass and Stiffness matrices
include("Assemble.jl")

# Time solver parameters
dt = CFL*dt
dtmin = dt
half_dt = 0.5*dtmin
half_dt_sq = 0.5*dtmin^2

dtmax = 50 * 24 * 60*60		# 50 days

# dt modified slightly for damping
if ETA != 0
	dt = dt/sqrt(1 + 2*ETA)
end

# Initialize kinematic field: global arrays
global d = zeros(nglob)
global v = zeros(nglob)
v[:] = 0.5e-3
global a = zeros(nglob)



#......................
# Boundary conditions : absorbing boundaries on 3 sides, fault boundary on one side
#......................

# Left boundary
impedance = rho1*vs1
BcLC, iBcL = BoundaryMatrix(wgll, NelX, NelY, iglob, 0.5*dye, impedance, 'L')

# Right Boundary = free surface: nothing to do

# Top Boundary
impedance = rho1*vs1
BcTC, iBcT = BoundaryMatrix(wgll, NelX, NelY, iglob, 0.5*dxe, impedance, 'T')

# Mass matrix at boundaries
Mq = M[:]
M[iBcL] .= M[iBcL] .+ half_dt*BcLC
M[iBcT] .= M[iBcT] .+ half_dt*BcTC


# Dynamic fault at bottom boundary
FltB, iFlt = BoundaryMatrix(wgll, NelX, NelY, iglob, 0.5*dxe, 1, 'B') # impedance = 1

const FltZ = M[iFlt]./FltB /half_dt * 0.5
const FltX = x[iFlt]


#.....................
# Initial Conditions 
#.....................
include("initialConditions/defaultInitialConditions.jl")
#tauo = repmat([22.5e6], FaultNglob, 1)

cca, ccb = fricDepth(cca, ccb, FltX)

Seff = SeffDepth(FltX)

tauo = tauDepth(FltX)

#.....................................
# Stresses and time related variables
#.....................................
tau = zeros(FaultNglob)

psi = tauo./(Seff.*ccb) - fo./ccb - (cca./ccb).*log.(2*v[iFlt]./Vo)
psi0 = psi[:]

Nel_ETA = 0
# Kelvin-Voigt Viscosity
if ETA !=0
	Nel_ETA = NelX
	x1 = 0.5*(1 + xgll')
	eta_taper = exp.(-pi*x1.^2)
	eta = ETA*dt * repmat([eta_taper], NGLL, 1)

else
    Nel_ETA = 0
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
tvsx = 2*yr2sec #	Output slip every 0.2 years
tvsxinc = tvsx
nevne = 0
tevneinc = 0.5
Vevne = Vthres

# Compute XiLf for each fault node
Xith = zeros(FaultNglob)
XiLf = zeros(FaultNglob)

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
OutLoc1 = 1e3
OutLoc2 = 2e3
OutLoc3 = 6e3
OutLoc4 = 8e3

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
println("Average node spacing = ", LX/(FaultNglob-1))
@printf("dt = %1.17f", dt)


# Find nodes that do not belong to fault
FltNI = deleteat!(collect(1:nglob), iFlt) 

# Some more initializations
r = zeros(nglob)
beta_ = zeros(nglob)
alpha_ = zeros(nglob)
p = zeros(nglob)

F = zeros(nglob)
dPre = zeros(nglob)
vPre = zeros(nglob)
dd = zeros(nglob)
dacum = zeros(nglob)

dnew = zeros(length(FltNI))


################################
# Time solver variables
################################

# Initially 1 = quasistatic; 2 = dynamic
isolver = 1

# Compute the diagonal of K
Kdiag = zeros(nglob)
Klocdiag = zeros(NGLL, NGLL)
for e = 1:Nel
    ig = iglob[:,:,e]
    wloc = W[:,:,e]
    Klocdiag[:,:] = 0

    for k =  1:NGLL
        for j = 1:NGLL
            Klocdiag[k,j] = Klocdiag[k,j] + 
            sum( coefint1[e]*H[k,:].*(wloc[:,j].*Ht[:,k])
                + coefint2[e]*(wloc[k,:].*H[j,:]).*Ht[:,j] )
        end
    end

    Kdiag[ig] .= Kdiag[ig] .+ Klocdiag[:,:]
end

diagKnew = Kdiag[FltNI]


v[:] = v[:] - 0.5*Vpl
Vf = 2*v[iFlt]
const iFBC = find(abs.(FltX) .> Fdepth)
const NFBC = length(iFBC)
Vf[iFBC] = 0

# Fault boundary: indices where fault within 24 km
fbc = reshape(iglob[:,1,:], length(iglob[:,1,:]),1)
idx = find(fbc .== find(x .>= -Fdepth)[1] - 1)[1]
const FltIglobBC = fbc[1:idx]

v[FltIglobBC] = 0

#idelevne = 0


# Output slip and recurrence interval at the start of every cycle
flag = 0
ev_it = 1
ev_it2 = 1
Coslip = zeros(FaultNglob, 1e3)
rec_int = zeros(1e3)


# Preallocate variables with unknown size
time_ = zeros(1e5)

delfsec = zeros(FaultNglob, 1e4)
Vfsec = zeros(FaultNglob, 1e4)
Tausec = zeros(FaultNglob, 1e4)

delf5yr = zeros(FaultNglob, 1e4)
Vf5yr = zeros(FaultNglob, 1e4)
Tau5yr = zeros(FaultNglob, 1e4)

Stress = zeros(FaultNglob, 1e5)
SlipVel = zeros(FaultNglob, 1e5)
Slip = zeros(FaultNglob, 1e5)
State = zeros(FaultNglob, 1e5)

println("\nSetup Complete")
