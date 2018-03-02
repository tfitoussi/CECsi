#/bin/python2

from numpy import loadtxt, size, pi, exp, logspace
from matplotlib.pyplot import figure, show, savefig
from analysisLib.analytic import *
from analysisLib.read import path_to_fig_dir

import matplotlib as mpl
label_size = 20
mpl.rcParams['xtick.labelsize'] = label_size 
mpl.rcParams['ytick.labelsize'] = label_size 
mpl.rcParams['ps.fonttype'] = 42

def nCMB(E,z):
   kTcmb = k*Tcmb*erg*(1+z)
   return (hb*c*erg)**(-3) *(E/pi)**2 /(exp(E/kTcmb)-1)    

redshifts       =[0,0.01,0.13,2]
# always take line number (column 0 for h.nu)
ind_dominguez   =[1,3,5,15]     # see z_Dominguez.dat
ind_lower_limit =[1,54,108,179] # see z_IR.dat 
ind_franceschini=[1,1,2,11]     # step of 0.2 (11 redshifts)
ind_finke       =[1,2,14,201]   # step of 0.01 (500 redshifts)
ind_gilmore     =[1,3,6,14]     # see z_IR_gil.dat 
ind_best_fit    =[1,1,6,81] 

i=0
for z in redshifts:
    ax = figure(figsize=(12,9)).add_subplot(111)
   
    #==== CMB ====
    hv = logspace(-4,-1,1000)
    ax.plot(hv,nCMB(hv,z)*hv**2,"-k",linewidth=2,label="CMB")

    # ==== Dominguez ====
    lamb,lambdaI = loadtxt("cascade_program/EBL_files/lambdaI_Dominguez.dat",unpack=True,usecols=[0,ind_dominguez[i]])
    hv = h*c/(lamb*1e-4)  # erg
    density = 1e-6*(4*pi/c)*lambdaI/(hv**2) /(erg) *(1+z)**3
    hv = hv *erg
    ax.plot(hv,density*hv**2,"-y",linewidth=2,label="Dominguez et Al 2011")
 
    # ==== Finke ====
    hv,density = loadtxt("cascade_program/EBL_files/n_Finke.dat",unpack=True,usecols=[0,ind_finke[i]])
    ax.plot(hv,density*(1+z)**3*hv**2,"-b",linewidth=2,label="Finke et Al 2010")
      
    # ==== Fransceschini ====
    hv,density = loadtxt("cascade_program/EBL_files/n_Fra_2008.dat",unpack=True,usecols=[0,ind_franceschini[i]])
    ax.plot(hv,density*(1+z)**3*hv**2,"-r",linewidth=2,label="Franceschini et Al 2008")
    hv,density = loadtxt("cascade_program/EBL_files/n_Fra_2017.dat",unpack=True,usecols=[0,ind_franceschini[i]])
    ax.plot(hv,density*(1+z)**3*hv**2,"--r",linewidth=2,label="Franceschini et Al 2017")

    # ==== Gilmore ====
    hv,density = loadtxt("cascade_program/EBL_files/n_Gil.dat",unpack=True,usecols=[0,ind_gilmore[i]])
    ax.plot(hv,density*hv**2,"-g",linewidth=2,label="Gilmore et Al 2012")
      
    # ==== Kneiske and Doll - "best fit" ====
    hv,density = loadtxt("cascade_program/EBL_files/n_bestfit10.dat",unpack=True,usecols=[0,ind_best_fit[i]])    
    ax.plot(hv,density*(1+z)**3*hv**2,"-c",linewidth=2,label="Kneiske et Al 2004 - 'best fit'")

    # ==== Kneiske and Doll - "lower limit" ====
    hv,density = loadtxt("cascade_program/EBL_files/n_lowerlimit10.dat",unpack=True,usecols=[0,ind_lower_limit[i]])
    ax.plot(hv,density*(1+z)**3*hv**2,"-m",linewidth=2,label="Kneiske et Al 2004 - 'lower limit'")

    ax.set_title("z=%.2f"%z,fontsize=20)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlim([1e-4,20])
    ax.set_ylim([1e-5,50])
    ax.grid(b=True,which='major')
    ax.legend(loc="best",fontsize=18)#,frameon=False,framealpha=0.5)
    ax.set_xlabel("energy [eV]",fontsize=20)
    ax.set_ylabel("$E^2 dN/dE$ [photon.eV.cm$^{-3}$]",fontsize=20)
    
    savefig(path_to_fig_dir+"ExtragalacticLight_z=%1.2f.eps"%z,bbox_inches='tight')  
      
    i+=1

    show()
