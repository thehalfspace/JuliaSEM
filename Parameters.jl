#######################################################################
#
#		PARAMETER FILE: SET THE PHYSICAL PARAMETERS FOR THE SIMULATION
#
#######################################################################


#----------------
#	Domain Size:
#----------------
distN = 1000	#	km to m conversion
Nsize = 2
LX = Nsize*24e3/distN	#	Length of Horizontal dimension of box	
LY = Nsize*15e3/distN	#	Length of Vertical dimension of box

NelX = 15	#	No. of elements in X
NelY = 10	#	No. of elements in Y

NelX = NelX*Nsize
NelY = NelY*Nsize

dxe = LX/NelX #	Size of one element along X
dye = LY/NelY #	Size of one element along Y

Nel = NelX*NelY # Total no. of elements


#----------------
#	No. of nodes
#----------------
P = 4		#	Lagrange polynomial degree
NGLL = P + 1 #	No. of Gauss-Legendre-Lobatto nodes

FaultNglob = NelX*(NGLL - 1) + 1

#---------------------------------
#	Parameters of the time solver
#---------------------------------
yr2sec = 365*24*60*60
Total_time = 20*yr2sec
CFL = 0.6	#	Courant stability number
dt = Inf	#	Timestep: set later


#------------------------------------
#	Jacobian for the global -> local 
#	coordinate conversion
#------------------------------------
dx_dxi = 0.5*dxe
dy_deta = 0.5*dye
jac = dx_dxi*dy_deta
coefint1 = jac/dx_dxi^2
coefint2 = jac/dy_deta^2

#-------------------------------------
#	Physical properties of the medium
#	(Modify them when assembling Mass
#	 and Stiffness matrices )
#-------------------------------------
rho1 = 2670
vs1 = 3464
rho2 = 2500
vs2 = 0.6*vs1
ETA = 0

rho = zeros(NGLL, NGLL)
mu = zeros(NGLL, NGLL)

# Low velocity layer dimensions
ThickX = 0
ThickY = 0


#--------------------------
# Earthquake parameters
#--------------------------
Vpl = 35e-3/yr2sec	#	Plate loading

Seff= repmat([50e6], FaultNglob, 1)		#	Effective normal stress
fo 	= repmat([0.6], FaultNglob, 1)		#	Reference friction coefficient
cca = repmat([0.015], FaultNglob, 1)	#	Rate-state parameter 'a'
ccb = repmat([0.020], FaultNglob, 1)	#	Rate-state parameter 'b'
Vo 	= repmat([1e-6], FaultNglob, 1)		#	Reference velocity 'Vo'
xLf = repmat([0.008/distN], FaultNglob, 1)#	Dc (Lc) = 8 mm
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
