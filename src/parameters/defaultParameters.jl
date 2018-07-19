#######################################################################
#
#		PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION
#
#######################################################################


#----------------
#	Domain Size:
#----------------
const distN = 1	#	km to m conversion
const Nsize = 2
const LX = Nsize*24e3/distN	#	Length of Horizontal dimension of box	
const LY = Nsize*15e3/distN	#	Length of Vertical dimension of box

const NelX = 15*Nsize	#	No. of elements in X
const NelY = 10*Nsize 	#	No. of elements in Y


# for constant mesh size using MeshBox, uncomment the two lines below
#dxe = LX/NelX #	Size of one element along X
#dye = LY/NelY #	Size of one element along Y

# for variable mesh size using VariableMesh, use the lines below
x_points = LX*sin.(pi/2*collect(0:NelX)/NelX)
y_points = LY*(1 - sin.(pi/2*collect(0:NelY)/NelY))
dxe = diff(x_points)
dye = abs.(diff(y_points))

const Nel = NelX*NelY # Total no. of elements

Fdepth = 24e3/distN # fault depth

#----------------
#	No. of nodes
#----------------
const P = 4		#	Lagrange polynomial degree
const NGLL = P + 1 #	No. of Gauss-Legendre-Lobatto nodes

const FaultNglob = NelX*(NGLL - 1) + 1

#---------------------------------
#	Parameters of the time solver
#---------------------------------
const yr2sec = 365*24*60*60
Total_time = 10*yr2sec
const CFL = 0.6	#	Courant stability number
dt = Inf	#	Timestep: set later

const IDstate = 2

#-------------------------------------
#	Physical properties of the medium
#	(Modify them when assembling Mass
#	 and Stiffness matrices )
#-------------------------------------
const rho1 = 2670
const vs1 = 3464
const rho2 = 2500
const vs2 = 0.6*vs1
ETA = 0

rho = zeros(NGLL, NGLL)
mu = zeros(NGLL, NGLL)

# Low velocity layer dimensions
const ThickX = 0e3
const ThickY = 0e3


#--------------------------
# Earthquake parameters
#--------------------------
const Vpl = 35e-3/yr2sec	#	Plate loading

Seff= repmat([50e6], FaultNglob)		#	Effective normal stress
tauo = repmat([22.5e6], FaultNglob)  #   Initial shear stress
fo 	= repmat([0.6], FaultNglob)		#	Reference friction coefficient
cca = repmat([0.015], FaultNglob)	#	Rate-state parameter 'a'
ccb = repmat([0.020], FaultNglob)	#	Rate-state parameter 'b'
const Vo 	= repmat([1e-6], FaultNglob)		#	Reference velocity 'Vo'
const xLf = repmat([0.008/distN], FaultNglob)#	Dc (Lc) = 8 mm
FaultC = zeros(FaultNglob)
Vf1 = zeros(FaultNglob)
Vf2 = zeros(FaultNglob)
Vf 	= zeros(FaultNglob)
psi1= zeros(FaultNglob)
psi2= zeros(FaultNglob)
tau1= zeros(FaultNglob)
tau2= zeros(FaultNglob)
tau3= zeros(FaultNglob)
tauNR= zeros(FaultNglob)
tauAB= zeros(FaultNglob)

#-----------------------
#	Output Seismograms
#-----------------------
OUTxseis = collect(-15:3:0)
OUTnseis = length(OUTxseis)
OUTyseis = repmat([15],OUTnseis)
