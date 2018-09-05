from numpy import sqrt, abs, cos, sin, arcsin, searchsorted, loadtxt, log, array
from numpy import select, exp, log10
from numpy import pi, sqrt, interp, dtype, zeros, size
from scipy.integrate import quad, fixed_quad
from CEClib.constants import *

# Value for analytic expression
Ee_default=1e3 #GeV
Egamma_default=1 #GeV
B_default=1e-15 #Gauss
lambda_B_default=0.3 #Mpc

def Analytic_dNdE(E_ic,zs=0):
   return  me*1e3/4*sqrt(3e9/Ecmb)*1e-9*E_ic**(-3/2.) / (1+zs) #GeV

def Analytic_dNdtheta(theta,B=B_default,tau_gg=397.4):
   return m*c**2 *sqrt((3*sigmaT*nCMB)/(2*e*Ecmb*eV) * tau_gg/(theta/degre*B)) # rad^-1

def Analytic_dNdt(t,B=B_default,lgg=1.32):
   return (m*c**2/4) *sqrt((3*sigmaT*nCMB)/(e*Ecmb*eV*B))*(8*c*t*yr/(lgg*Mpc))**(0.25) /t

def Analytic_observables_vs_energy(Egamma,zs,B=B_default,Eg0=1e5):
   Ds = distance(zs)[1]
   Ee0=Eg0*1e3/2 #GeV
   RL0=RL(Ee0,B,zs)
   Dic0=Dic(Ee0,zs)
   E_e = Ee(Egamma,zs)
   delta = Dic0/(2*RL0)*((Ee0/E_e)**2 -1)
   lgg = 1.4 #sslambda_gg(Eg0*1e3,zs)[1]
   theta = arcsin(lgg/Ds*sin(delta))
   c_delta_t = lgg*(1-cos(delta)) - Ds*(1-cos(theta)) #+Dic(E_e,zs)
   lgg = 117 #lambda_gg(8e3,zs)[1]
   theta2 = arcsin(lgg/Ds*sin(delta))
   c_delta_t2 = lgg*(1-cos(delta)) - Ds*(1-cos(theta)) #+Dic(E_e,zs)
   return abs(theta)/degre, c_delta_t *Mpc/c, abs(theta2)/degre, c_delta_t2 *Mpc/c # sec.

def Analytic_delay_vs_theta(theta,zs,Esource=1e5):
   distSource=distance(zs)[1]
   E0_source = Esource*1e-3 #TeV
   Dic0=Dic(E0_source/2)
   lgg = 1.32#lambda_gg(Esource*1e3)[0]
   delta_ic = arcsin(distSource/lgg*sin(theta))
   c_delta_t = (lgg*(1-cos(delta_ic)) - distSource*(1-cos(theta)))
   lgg = 117#lambda_gg(8e3)[0]
   delta_ic = arcsin(distSource/lgg*sin(theta))
   c_delta_t2 = (lgg*(1-cos(delta_ic)) - distSource*(1-cos(theta)))
   return c_delta_t *Mpc/c, c_delta_t2 *Mpc/c  # sec.

def arrival_approx(B=B_default,z=0.13,Egamma0=1e5): # Pour Egamma = 1GeV
   delta = Dic(Ee=Ee(1),z=z)/(2*RL(B=B,z=z))  
   lgg =  1.32# lambda_gg(Egamma0,z)[0]
   #print distance(z)[0], lgg, distance(z)[0]/lgg
   theta = lgg/distance(z)[0] * delta
   cdt = lgg/2 * delta**2
   return theta/degre, cdt*Mpc/c /yr

# Compton accumulation
def ECompton_threshold(Compton_threshold = 0.005,z=0):
   return Compton_threshold/(4/3*Ecmb*(1+z)/me*1e-3) *me*1e-6 #GeV

# Compton scattering
def Dic(Ee=Ee_default,z=0): # Ee (GeV)
   return 3*(me*1e3)**2/(4*sigmaT*Ee*1e9*rhoCMB*(1+z)**4) /Mpc #Mpc

def lambdaIC(z=0):
   return 1/(sigmaT*Mpc*nCMB*(1+z)**3) #Mpc

def Eic(Ee=Ee_default,z=0):
   return 4*Ecmb*(1+z) *Ee**2/(3*me**2)*1e3 #GeV 

def Ee(Egamma=Egamma_default,z=0):
   return me*sqrt((3*Egamma*1e-3 )/(4*Ecmb*(1+z))) #GeV 

def tIC():
   return lambdaIC()/(c*yr/Mpc) #yr

# Larmor radius
def RL(Ee=Ee_default,B=B_default,z=0):
   return (Ee*GeV)/(e*B) /(1+z) /Mpc #Mpc

# Magnetic deflection
def delta(Ee=Ee_default,B=B_default,lambda_B=lambda_B_default,z=0):
   Ee=array(Ee)
   D_ic=Dic(Ee,z)
   condlist=[lambda_B>=D_ic,lambda_B<D_ic]
   delta=[D_ic/RL(Ee,B),sqrt(D_ic*lambda_B)/RL(Ee,B)]
   return select(condlist,delta)

def Delta(Ee=Ee_default,B=B_default,lambda_B=lambda_B_default,z=0,Esource=1e5):
   lgg=lambda_gg(Esource)[0]
   Dsource=distance(z)[1]
   return arcsin(lgg/Dsource*delta(Ee,B,lambda_B,z))/degre

def Dt(Ee=Ee_default,B=B_default,lambda_B=lambda_B_default,z=0,Esource=1e5):
   lgg=lambda_gg(Esource)[0]
   return lgg/2*delta(Ee,B,lambda_B,z)**2 /c *Mpc/yr

#synchrotron radiation
def Esc(Ee=Ee_default,B=B_default): # eV
   return (3*hb*e*c)/(2*m**3*c**6) * B*(Ee*GeV)**2 *erg*1e9

# Threshold energy when Dic = RL
def Ethreshold_ic(Ee=Ee_default,B=B_default):
   return Eic(Ee)*Dic(Ee)/RL(Ee,B) # GeV

def Ethreshold_gg(epsilon=Eebl):
   return (me)**2/epsilon *1e-3 #GeV

def lambda_gg(Egamma=1e3,z=0): # Egamma (GeV)
   # Bilinear interpolation 
   z_tab = [0,0.03,0.1,0.5,1,2,3]
   E_tab = loadtxt("CEClib/lambda_e.dat",unpack=True,usecols=[0])*me*1e-6
   i2 = searchsorted(z_tab,z)
   if z<=0:
      fy=1
      i1=i2
   elif z>3:
      fy=1
      i2-=1
      i1=i2
   else:
      i1 = i2-1
      fy=(z-z_tab[i1])/(z_tab[i2]-z_tab[i1])
   j2 = searchsorted(E_tab,Egamma)
   j1 = j2-1
   fx=(Egamma-E_tab[j1])/(E_tab[j2]-E_tab[j1])
   lambda_e = loadtxt("CEClib/lambda_e.dat",unpack=True,usecols=[i1+1,i2+1,i1+7,i2+7])
   lambda11 = lambda_e[0,j1]*((1-fx)*(1-fy))
   lambda12 = lambda_e[1,j1]*((1-fx)*fy)
   lambda21 = lambda_e[0,j2]*(fx*(1-fy))
   lambda22 = lambda_e[1,j2]*(fx*fy)
   lambda_proper = lambda11 + lambda12 + lambda21 + lambda22
   lambda11 = lambda_e[2,j1]*((1-fx)*(1-fy))
   lambda12 = lambda_e[3,j1]*((1-fx)*fy)
   lambda21 = lambda_e[2,j2]*(fx*(1-fy))
   lambda22 = lambda_e[3,j2]*(fx*fy)
   lambda_comobile = lambda11 + lambda12 + lambda21 + lambda22
   return lambda_proper, lambda_comobile, 800.e3/Egamma #Mpc (from Durrer and Neronov 2013)

def comobileTime(z):
   return -1/(H0*(1+z)*sqrt(omegaM*(1+z)**3+omegaK*(1+z)**2+omegaL))

def distPhoton(z):
   return -c/(H0*a0*sqrt(omegaM*(1+z)**3+omegaK*(1+z)**2+omegaL))

def distLepton(z,E):
   beta = sqrt(1-m**2*c**4/E**2)
   return -beta*c/(H0*a0*sqrt(omegaM*(1+z)**3+omegaK*(1+z)**2+omegaL))

def properIntegrand(z):

   return -c/(H0*(1+z)*sqrt(omegaM*(1+z)**3+omegaK*(1+z)**2+omegaL))

def comobileIntegrand(z):
   return -c/(H0*sqrt(omegaM*(1+z)**3+omegaK*(1+z)**2+omegaL))

def distance(z=0.,Nquad=30):
   if type(z) == dtype(float): 
      distpr =  fixed_quad(properIntegrand,z,0,n=Nquad)[0]/Mpc
      distco =  fixed_quad(comobileIntegrand,z,0,n=Nquad)[0]/Mpc
   else:
      distpr = zeros(size(z))
      distco = zeros(size(z))
      for i in range(size(z)): 
         distpr[i] =  fixed_quad(properIntegrand,z[i],0,n=Nquad)[0]/Mpc
         distco[i] =  fixed_quad(comobileIntegrand,z[i],0,n=Nquad)[0]/Mpc
   return distco, distpr

def Ecut(z=0.14): # TeV
   zi,ctgg,best_fit,dominguez,finke,franceschini,gilmore,lower_limit = loadtxt("CEClib/Ecut.dat",unpack=True,usecols=[0,1,2,3,4,5,6,7]) 
   return interp(z,zi,best_fit),interp(z,zi,dominguez),interp(z,zi,finke),interp(z,zi,franceschini),interp(z,zi,gilmore),interp(z,zi,lower_limit)

def Etarget(Egamma=1): # Egamma en TeV, Etarget en eV
   return 6.9*(Egamma/0.2)**(-0.9)

def PSF_Taylor2011(E): # E in GeV -> PSF in degre
   return 1.7 * E**(-0.74) * (1+(E/15)**2)**0.37

def intrinsic_luminosity(Gamma,Emin,Emax):
   if Gamma == 1:
      return (Emax**2-Emin**2)/(Emin*log(Emax/Emin))
   elif Gamma == 2:
      return log(Emax/Emin)/(1/Emin-1/Emax)
   else:
      return (1-Gamma)/(2-Gamma) * (Emax**(2-Gamma) - Emin**(2-Gamma))/(Emax**(1-Gamma) - Emin**(1-Gamma))
    
def norma_spectrum(Gamma,Emin,Emax):
   if Gamma == 1:
      return 1/log(Emax/Emin)
   else:
      return (1-Gamma)/((Emax/Emin)**(1-Gamma)-1) 
