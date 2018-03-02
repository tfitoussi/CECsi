from scipy.special import zeta
from numpy import pi, sqrt

e=4.8032068e-10 # esu (cgs units)
qe=1.602176565e-12 # erg
m=9.1093897e-28 # g
me=510.998928 # kev
c=2.99792458e10 # cm.s-1
k=1.380658E-16 # erg.K-1
h=6.62606957e-27 # erg.s
hb=1.05457266E-27 # erg.s
alpha=1./137.035999679
lambdaC=hb/(m*c)

# conversions
Mpc=(3.0856776e+16)*1e8 # Mpc to cm
eV=1.602176565e-12 # erg
keV = 1e3  *eV # erg
MeV = 1e6  *eV # erg
GeV = 1e9  *eV # erg
TeV = 1e12 *eV # erg
erg = 1/eV # eV
degre = pi/180 # rad 
rad = 180/pi # degre 
day=24*3600 # s
yr=365.25*day # s

# Cosmology
a0=1.
hubble = 0.678 # Hubble constant
H0=100*(hubble*1e5)/(Mpc) # s-1
omegaM = 0.3
omegaK = 0
omegaL = 0.7

# CMB 
Tcmb=2.7255 # K
Tcmbp=(k*Tcmb)/(m*c*c) # Adimentionnal thermal energy
nCMB=(16*pi*zeta(3,1))*(k*Tcmb/(h*c))**3 # cm^-3
Ecmb=(pi**4/(30*zeta(3,1)))*Tcmbp*me*1e3 #eV
rhoCMB=Ecmb*nCMB #eV.cm^-3

# EBL
Eebl=1 #eV

# Particles physics
r0=e*e/(m*c*c) # cm
sigmaT=8.*pi*(r0**2)/3. # cm^2, Thomson cross section
