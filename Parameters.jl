#######################################################################
#
#		PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION
#
#######################################################################


#----------------
#	Domain Size:
#----------------
const distN = 1000	#	km to m conversion
const Nsize = 2
const LX = Nsize*24e3/distN	#	Length of Horizontal dimension of box	
const LY = Nsize*15e3/distN	#	Length of Vertical dimension of box

const NelX = 15*Nsize	#	No. of elements in X
const NelY = 10*Nsize 	#	No. of elements in Y

#NelX = NelX*Nsize
#NelY = NelY*Nsize

const dxe = LX/NelX #	Size of one element along X
const dye = LY/NelY #	Size of one element along Y

const Nel = NelX*NelY # Total no. of elements


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
const Total_time = 0.1*yr2sec
const CFL = 0.6	#	Courant stability number
dt = Inf	#	Timestep: set later


#------------------------------------
#	Jacobian for the global -> local 
#	coordinate conversion
#------------------------------------
const dx_dxi = 0.5*dxe
const dy_deta = 0.5*dye
const jac = dx_dxi*dy_deta
const coefint1 = jac/dx_dxi^2
const coefint2 = jac/dy_deta^2

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
const ThickX = 0
const ThickY = 0


#--------------------------
# Earthquake parameters
#--------------------------
const Vpl = 35e-3/yr2sec	#	Plate loading

Seff= repmat([50e6], FaultNglob, 1)		#	Effective normal stress
fo 	= repmat([0.6], FaultNglob, 1)		#	Reference friction coefficient
cca = repmat([0.015], FaultNglob, 1)	#	Rate-state parameter 'a'
ccb = repmat([0.020], FaultNglob, 1)	#	Rate-state parameter 'b'
const Vo 	= repmat([1e-6], FaultNglob, 1)		#	Reference velocity 'Vo'
const xLf = repmat([0.008/distN], FaultNglob, 1)#	Dc (Lc) = 8 mm
FaultC = zeros(FaultNglob, 1)
Vf1 = zeros(FaultNglob, 1)
Vf2 = zeros(FaultNglob, 1)
Vf 	= zeros(FaultNglob, 1)
psi1= zeros(FaultNglob, 1)
psi2= zeros(FaultNglob, 1)
tau1= zeros(FaultNglob, 1)
tau2= zeros(FaultNglob, 1)
tau3= zeros(FaultNglob, 1)
tauNR= zeros(FaultNglob, 1)
tauAB= zeros(FaultNglob, 1)

#-----------------------
#	Output Seismograms
#-----------------------
OUTxseis = collect(-15:3:0)
OUTnseis = length(OUTxseis)
OUTyseis = repmat([15],OUTnseis,1)
