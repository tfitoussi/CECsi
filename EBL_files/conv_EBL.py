#!/bin/python

from numpy import loadtxt, arange, zeros, savetxt, interp, log10, size, shape, pi
from sys import argv

eV=1.602176565e-12 # erg
erg = 1/eV # eV
c=2.99792458e10 # cm.s-1
h=6.62606957e-27 # erg.s

if argv[1] == "Dominguez":
    data = loadtxt("Dominguez_2011_lambdaI.dat")
    lamb = data[:,0]
    lambdaI = data[:,1:]
    z = loadtxt("Dominguez_2011_z.dat")
    hv = h*c/(lamb*1e-4)  # erg
    
    density = zeros((size(hv)+1,size(z)+1))
    density[1:,0] = hv[::-1] *erg
    density[0,1:] = z
    for i in range(size(z)-1):
        density[1:,size(z)-(i+1)] = 1e-6*(4*pi/c)*lambdaI[:,i]/(hv**2) /(erg) *(1+z[i])**3
    
    savetxt("Dominguez_2011.dat",density[:,:],fmt='%1.4e')

if argv[1] == "Franceschini":
    years = ["2008","2017"]
    for year in years:
        data = loadtxt("Franceschini_"+year+"_n.dat")
    
        E_eV = 10**data[:,0]
    
        density = zeros((32,12))
        dens = zeros(31)
        density[1:,0] = E_eV
        j=1
        z=0.2*arange(11)
        density[0,1:] = z
        for i in arange(1,22,2):
            dens = 10**(data[:,i])/(E_eV*(1+z[j-1]))
            density[1:,j] = interp(log10(E_eV),log10(E_eV*(1+z[j-1])),dens) / (1+z[j-1])**3
            j+=1
    
        savetxt("Franceschini_"+year+".dat",density[:,:],fmt='%1.4e')


