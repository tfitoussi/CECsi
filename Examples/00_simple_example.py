#!/usr/bin/python3
#-*- coding: utf-8 -*-

#=================================== 00 Simple Example ===================================#
# This simple example aims to show how-to simple launch a simulation of cosmological
# electromagnetic cascade and extract data 
#=========================================================================================#

from CECsimulation import *
from distribution_graph import *

def simple_simulation(redshift,EGMF_B):
    simu = CECsimulation(redshift,EGMF_B)
    #fig1 = distribution_graph()
    fig2 = distribution_graph()
    #fig3 = distribution_graph()

    # all generations
    simu.Run()
    #E, dNdE = simu.Distribution("energy")
    #fig1.Add_histo(E,E**2 * dNdE * (1+redshift),"total",color="k")
    theta, dNdtheta = simu.Distribution("angle")
    fig2.Add_histo(theta,theta *dNdtheta,"total",color="k")
    #dt, dNdt = simu.Distribution("time")
    #fig3.Add_histo(dt,dt *dNdt,"total",color="k")

    # plot generation by generation
    color=["b","r","g"]
    for gen in [1,2,3]:
        simu.Run()
        simu.Data_selection(generation=gen)
        #E, dNdE = simu.Distribution("energy")
        #fig1.Add_histo(E,E**2 * dNdE * (1+redshift),"gen = %d"%gen,color=color[gen-1])
        theta, dNdtheta = simu.Distribution("angle")
        fig2.Add_histo(theta,theta *dNdtheta,"gen = %d"%gen,color=color[gen-1])
        #dt, dNdt = simu.Distribution("time")
        #fig3.Add_histo(dt,dt *dNdt,"gen = %d"%gen,color=color[gen-1])

    #fig1.Plot(x_range=[1e-3,1e5],y_range=[1e2,5e4],xlabel="$E$ [GeV]",ylabel="$E^2 dN/dE$ [GeV]",leg1_pos="upper right")
    fig2.Plot(x_range=[1e-5,180],y_range=[1e0,2e3],xlabel="$\\theta$ [deg.]",ylabel="$\\theta dN/d\\theta$ [/ primary phot.]",leg1_pos="lower left")
    #fig3.Plot(x_range=[1e-5,1e10],y_range=[1e-1,2e3],xlabel="$\\delta t$ [yr]",ylabel="$\\delta t dN/dt$ [/ primary phot.]",leg1_pos="lower left")



if __name__ == '__main__':
    simple_simulation(redshift=0.13,EGMF_B=3e-16)
