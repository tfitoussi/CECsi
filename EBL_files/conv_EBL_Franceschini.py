#!/bin/python

from numpy import loadtxt, arange, zeros, savetxt, interp, log10

years = ["2008","2017"]

for year in years:
    data = loadtxt("Fra_"+year+".dat")

    E_eV = 10**data[:,0]

    density = zeros((31,12))
    dens = zeros(31)
    density[:,0] = E_eV
    j=1
    z=0.2*arange(11)
    for i in arange(1,22,2):
        dens = 10**(data[:,i])/(E_eV*(1+z[j-1]))
        density[:,j] = interp(log10(E_eV),log10(E_eV*(1+z[j-1])),dens) / (1+z[j-1])**3
        j+=1

    savetxt("n_Fra_"+year+".dat",density[::-1,:],fmt='%1.4e')
