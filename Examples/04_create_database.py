#!/usr/bin/python3
#-*- coding: utf-8 -*-

#============================ 04 Create a simulation database ============================#
# Analysing multiple different cases (of EGMF for example) requires to run multiple MC
# simulations. This can take a long time. Thus, it is useful to generate a database of
# this simulations first. 
# Here we generate the simulation for different redshifts, amplitudes and coherence length
# of the EGMF with an intrinsic spectrum of the source extenting from 1 GeV to 1 TeV.
#=========================================================================================#

from CECsimulation import *
from numpy import logspace

def Generate_DB(redshift_list=[0.14],EGMF_B_list=[1e-15],EGMF_LB_list=[1]):
    for z in redshift_list:
        for B in EGMF_B_list:
            for LB in EGMF_LB_list:
                simu = CECsimulation(z,B,LB,Emin=1e-3,Emax=100)
                simu.Run()


if __name__ == '__main__':
    redshift_list = [0.04,0.14,0.4,1,2]
    B_list = logspace(-22,-12,11)
    LB_list = logspace(-5,3,9)
    Generate_DB(redshift_list,B_list,LB_list)
