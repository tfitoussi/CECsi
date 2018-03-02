#!/usr/bin/python3
#-*- coding: utf-8 -*-
from numpy import append
from matplotlib.pyplot import figure, show, savefig, legend, setp, gca, sca

from matplotlib import rcParams


label_size = 30
rcParams['xtick.labelsize'] = label_size
rcParams['ytick.labelsize'] = label_size
rcParams['ps.fonttype'] = 42

class GraphDistribution(object):
    def __init__(self,xlabel="",ylabel="",title=""):
        self.ax = figure(figsize=(12,9),tight_layout=True).add_subplot(111)
        self.ax.set_xscale('log')
        self.ax.set_yscale('log')
        self.ax.grid(b=True,which='major')

        self.ax.set_xlabel(xlabel,fontsize=label_size)
        self.ax.set_ylabel(ylabel,fontsize=label_size)
        if not (title == ""):
            self.ax.set_title(title,fontsize=label_size)
 
        self.legend1 = []
        self.plotlist1 = []
        self.legend2 = []
        self.plotlist2 = []

    def Add(self,x0,y0,label="",leg=1,linestyle="-",linewidth="2",color="k",histo=True):
        if histo:
            plot = self.ax.plot(x0,y0,color=color,linewidth=linewidth,linestyle=linestyle,drawstyle='steps-mid')
        else:
            plot = self.ax.plot(x0,y0,color=color,linewidth=linewidth,linestyle=linestyle)
        if leg==1:
            self.legend1 = append(self.legend1,[label])
            self.plotlist1 = append(self.plotlist1,[plot])
        else:
            self.legend2 = append(self.legend2,[label])
            self.plotlist2 = append(self.plotlist2,[plot])

    def Add_data(self,Xmed,Xmax,Xmin,Ymed,Ymax,Ymin,label,leg=1,marker='o',color="k",linewidth="2"):
        Xerr = [Xmed-Xmin,Xmax-Xmed]
        Yerr = [Ymax-Ymed,-Ymin+Ymed]
        plot = self.ax.errorbar(Xmed,Ymed,xerr=Xerr,yerr=Yerr,linestyle="",linewidth=linestyle,color=color,marker=marker)
        if leg==1:
            self.legend1 = append(self.legend1,[label])
            self.plotlist1 = append(self.plotlist1,[plot])
        else:
            self.legend2 = append(self.legend2,[label])
            self.plotlist2 = append(self.plotlist2,[plot])

    def Plot(self,x_range=[],y_range=[],leg1_pos="best",leg2_pos="best",save=""):
        sca(self.ax)
        if self.legend1 != []:
            leg = legend(self.plotlist1,self.legend1,fontsize='xx-large',loc=leg1_pos)
            gca().add_artist(leg)
        if self.legend2 != []:
            leg = legend(self.plotlist2,self.legend2,fontsize='xx-large',loc=leg2_pos)
            gca().add_artist(leg)
        
        if x_range != []:
            self.ax.set_xlim(x_range)
        if y_range != []:
            self.ax.set_ylim(y_range)

        if not (save == ""):
            savefig(save,bbox_inches='tight')
