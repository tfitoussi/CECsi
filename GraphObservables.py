#!/usr/bin/python3
#-*- coding: utf-8 -*-
from numpy import histogram2d, log10, flipud, rot90, ma, logspace, zeros, arange
from matplotlib.pyplot import figure, show, savefig, legend, setp, gca
from matplotlib.gridspec import GridSpec
from matplotlib.colorbar import Colorbar
import matplotlib.cm as colormap
from matplotlib import rcParams

label_size = 30
rcParams['xtick.labelsize'] = label_size
rcParams['ytick.labelsize'] = label_size
rcParams['ps.fonttype'] = 42

class GraphObservables(object):
    def __init__(self,title="",one_figure=True):
        self.one_figure = one_figure
        if one_figure:
            self.fig = figure(figsize=(20,18))
            gs = GridSpec(2, 2, height_ratios=[1,1], width_ratios=[1,1]) 
            self.fig.subplots_adjust(hspace=0,wspace=0)
            self.ax0 = self.fig.add_subplot(gs[1])
            self.ax1 = self.fig.add_subplot(gs[0])
            self.ax2 = self.fig.add_subplot(gs[2],sharex=self.ax1)
            self.ax3 = self.fig.add_subplot(gs[3],sharey=self.ax2)
        else:
            self.fig1 = figure(figsize=(12,9),tight_layout=True)
            self.ax1  = self.fig1.add_subplot(111)
            self.fig2 = figure(figsize=(12,9),tight_layout=True)
            self.ax2  = self.fig2.add_subplot(111)
            self.fig3 = figure(figsize=(12,9),tight_layout=True)
            self.ax3  = self.fig3.add_subplot(111)

        self.ax1.set_title(title,fontsize=label_size)
        self.ax1.yaxis.set_ticks([1e-7,1e-5,1e-3,1e-1,1e1,1e3,1e5,1e7,1e9])
        self.ax1.set_xscale('log')
        self.ax1.set_yscale('log')
        self.ax1.set_ylabel("$\\Delta t$ [yrs]",fontsize=label_size)
        if not (title == ""):
            self.ax1.set_title(title,fontsize=label_size)
        if one_figure:
            setp(self.ax1.get_xticklabels(), visible=False)
            self.ax0.axis('off')
        else:
            self.ax1.set_xlabel("$E$ [GeV]",fontsize=label_size)

        self.ax2.set_xscale('log')
        self.ax2.set_yscale('log')
        self.ax2.set_xlabel("$E$ [GeV]",fontsize=label_size)
        self.ax2.set_ylabel("$\\theta$ [deg]",fontsize=label_size)
 
        self.ax3.xaxis.set_ticks([1e-7,1e-5,1e-3,1e-1,1e1,1e3,1e5,1e7,1e9])
        self.ax3.set_xscale('log')
        self.ax3.set_yscale('log')
        self.ax3.set_xlabel("$\\Delta t$ [yrs]",fontsize=label_size)
        if one_figure:
            setp(self.ax3.get_yticklabels(), visible=False)
        else:
            self.ax3.set_ylabel("$\\theta$ [deg]",fontsize=label_size)

    
    def AddDensity(self,energy,theta,delay,weight,nbBins=250,cmap=colormap.YlOrBr):
        cond = (theta>0) & (delay>0)
        H, xedges, yedges = histogram2d(log10(energy[cond]),log10(delay[cond]),bins=nbBins,weights=weight[cond])
        H = flipud(rot90(ma.masked_where(H==0,H)))
        self.im1 = self.ax1.pcolormesh(10**xedges,10**yedges,log10(H),cmap=cmap)   

        H, xedges, yedges = histogram2d(log10(energy[cond]),log10(theta[cond]),bins=nbBins,weights=weight[cond])
        H = flipud(rot90(ma.masked_where(H==0,H)))
        self.im2 = self.ax2.pcolormesh(10**xedges,10**yedges,log10(H),cmap=cmap)      
    
        H, xedges, yedges = histogram2d(log10(delay[cond]),log10(theta[cond]),bins=nbBins,weights=weight[cond])
        H = flipud(rot90(ma.masked_where(H==0,H)))
        self.im3 = self.ax3.pcolormesh(10**xedges,10**yedges,log10(H),cmap=cmap)


    def AddMean(self,energy,theta,delay,weight,nbBins=28,color="k",label=""):
        ener,angle,dt = self.observables_vs_energy(energy,theta,delay,weight,nbBins=nbBins,median=False)
        self.ax1.plot(ener,dt,color=color,linewidth=2,label=label)
        self.ax2.plot(ener,angle,color=color,linewidth=2,label=label)
        dt,angle,ener = self.observables_vs_delay(energy,theta,delay,weight,nbBins=nbBins,median=False)
        self.ax3.plot(dt,angle,color=color,linewidth=2,label=label)


    def AddCurve_dt_vs_E(self,energy,delay,color="k",linestyle="--",label=""):
        self.ax1.plot(energy,delay,color=color,linewidth=2,linestyle=linestyle,label=label)

    def AddCurve_theta_vs_E(self,energy,theta,color="k",linestyle="--",label=""):
        self.ax2.plot(energy,theta,color=color,linewidth=2,linestyle=linestyle,label=label)

    def AddCurve_theta_vs_dt(self,delay,theta,color="k",linestyle="--",label=""):
        self.ax3.plot(delay,theta,color=color,linewidth=2,linestyle=linestyle,label=label)


    def Plot(self,energy_range=[],theta_range=[],delay_range=[],save="",colorbar=True):
        self.ax1.grid(b=True,which='major')
        self.ax2.grid(b=True,which='major')
        self.ax3.grid(b=True,which='major')

        if energy_range != []:
            self.ax1.set_xlim(energy_range)
            self.ax2.set_xlim(energy_range)
        if theta_range != []:
            self.ax2.set_ylim(theta_range)
            self.ax3.set_ylim(theta_range)
        if delay_range != []:
            self.ax1.set_ylim(delay_range)
            self.ax3.set_xlim(delay_range)

        if self.one_figure:
            if colorbar:
                self.cbar1=self.fig.colorbar(self.im1, ax=self.ax0)
                self.cbar1.ax.set_ylabel("counts [log]",fontsize=label_size)

            self.ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,fontsize=label_size)

        else:
            if colorbar:
                self.cbar1=self.fig1.colorbar(self.im1, ax=self.ax1)
                self.cbar1.ax.set_ylabel("counts [log]",fontsize=label_size)
                self.cbar2=self.fig2.colorbar(self.im2, ax=self.ax2)
                self.cbar2.ax.set_ylabel("counts [log]",fontsize=label_size)
                self.cbar3=self.fig3.colorbar(self.im3, ax=self.ax3)
                self.cbar3.ax.set_ylabel("counts [log]",fontsize=label_size)

            self.ax1.legend(loc="best",fontsize=label_size)
            self.ax2.legend(loc="best",fontsize=label_size)
            self.ax3.legend(loc="best",fontsize=label_size)

        if save != "":
            savefig(save,bbox_inches='tight')




    def observables_vs_energy(self,energy,theta,delay,weight,Erange=[1e-3,1e5],median=False,nbBins=-1):#28):
         # nbBins=16, Erange=[1e-1,1e3] <=> Arlen 2014 (energy binning)
         if (nbBins==-1):
             bins=append(logspace(log10(Erange[0]),log10(0.9),100),logspace(0,log10(Erange[1]),16))
             nbBins=size(bins)-1
         if (nbBins==-2): # <=> Arlen 2014 (energy binning)
             nbBins=16 
             Erange=[1e-1,1e3]
             bins=logspace(log10(Erange[0]),log10(Erange[1]),nbBins+1)
         else:
             bins=logspace(log10(Erange[0]),log10(Erange[1]),nbBins+1)
         ener = (bins[1:nbBins+1]*bins[0:nbBins])**0.5
         #dtheta, dE = histogram(energy,bins,weights=weight*theta)
         #dt,     dE = histogram(energy,bins,weights=weight*delay)
         #dN,     dE = histogram(energy,bins,weights=weight)
         dt = zeros(nbBins,dtype="float64")
         dtheta = zeros(nbBins,dtype="float64")
         if median==False:
             for i in arange(nbBins):
                mask = (energy >= bins[i]) & (energy < bins[i+1])
                w = sum(weight[mask])
                t = sum(weight[mask]*delay[mask])
                angle = sum(weight[mask]*theta[mask])
                if w !=0:
                    dt[i]=t/w
                    dtheta[i]=angle/w
         else:
             for i in arange(nbBins):
                 mask = (energy >= bins[i]) & (energy < bins[i+1])
                 if True in mask:
                     dt[i] = weighted_median(delay[mask],weight[mask])
                     dtheta[i] = weighted_median(theta[mask],weight[mask])
                 else:
                     dt[i] = 0
                     dtheta[i] = 0
         return ener, dtheta, dt
    
    def observables_vs_theta(self,energy,theta,delay,weight,theta_range=[1e-5,1e2],median=False,nbBins=-1):#28):
        # <=> Taylor 2011 (time binning)
        cond = energy >0.1
        energy =energy[cond]
        theta=theta[cond]
        delay=delay[cond]
        weight=weight[cond]
        bins=logspace(log10(theta_range[0]),log10(theta_range[1]),nbBins+1)
        th = (bins[1:nbBins+1]*bins[0:nbBins])**0.5
        dE = zeros(nbBins,dtype="float64")
        dt = zeros(nbBins,dtype="float64")
        for i in arange(nbBins):
            mask = (theta >= bins[i]) & (theta < bins[i+1])
            w = sum(weight[mask])
            E = sum(weight[mask]*energy[mask])
            time = sum(weight[mask]*theta[mask])
            if w !=0:
                dE[i]=E/w
                dt[i]=time/w
        return dE, th, dt

    def observables_vs_delay(self,energy,theta,delay,weight,time_range=[10**-4,10**8.5],median=False,nbBins=100):#9): 
        # <=> Taylor 2011 (time binning)
        bins=logspace(log10(time_range[0]),log10(time_range[1]),nbBins+1)
        time = (bins[1:nbBins+1]*bins[0:nbBins])**0.5
        dE = zeros(nbBins,dtype="float64")
        dtheta = zeros(nbBins,dtype="float64")
        if median==False:
            for i in arange(nbBins):
                mask = (delay >= bins[i]) & (delay < bins[i+1])
                w = sum(weight[mask])
                E = sum(weight[mask]*energy[mask])
                angle = sum(weight[mask]*theta[mask])
                if w !=0:
                    dE[i]=E/w
                    dtheta[i]=angle/w
        else:
            for i in arange(nbBins):
                mask = (delay >= bins[i]) & (delay < bins[i+1])
                if True in mask:
                    dE[i] = weighted_median(energy[mask],weight[mask])
                    dtheta[i] = weighted_median(theta[mask],weight[mask])
                else:
                    dE[i] = 0
                    dtheta[i] = 0
        return time, dtheta, dE


