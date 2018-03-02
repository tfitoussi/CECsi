#!/usr/bin/python3
#-*- coding: utf-8 -*-
from os import chdir, getcwd, path, makedirs
import subprocess, shutil
from time import ctime
from numpy import loadtxt, log10, interp, append, histogram
from CEClib.make_path_dir import make_path_dir
from CEClib.write_parameters import write_preprocessing_f95, write_input_parameters_f95
from CEClib.constants import yr, rad, degre


# Define usefull path
pathname = path.dirname(__file__)        
CECsi_PATH = path.abspath(pathname)
resultsFile="/results.dat"
resultsFile_z_0="/results_z_0.dat"

class CECsimulation(object):
    def __init__(self,redshift,EGMF_B,EGMF_LB=1,Emin=100,Emax=100,Emin_simu=1e-3,EBLmodel="Dominguez",simu_path=""):
        '''
        Initialize parameters for a cosmological electromagnetic cascade.
        Mandatory: 
        * redshift -> source redshift
        * EGMF_B [G] -> amplitude of the extragalatic magnetic field
        Optional:
        * EGMF_LB [Mpc] -> coherence length of the extragalactic magnetic field
        * Emin & Emax [TeV] -> energy min and max of the source photons
        * Emin_simu [TeV] -> minimum energy of the photons recorded
        * EBLmodel -> EBL model can be choosed among
            "Dominguez" from Dominguez et al 2011
            "Franceschini" from Franceschini et al 2017
            "Finke" from Finke et al 2010
            "Gilmore" from Gilmore et al 2012
            "Lower limit" or "Best fit" from Kneiske et Al 2004   
        * simu_path allows to give the path to the simulation files (be careful when
        using it)
        '''
        # source parameters 
        self.redshift = redshift
        self.Emin = Emin # TeV 
        self.Emax = Emax # TeV 
        # Intergalactic medium parameters
        self.EBL = EBLmodel
        self.EGMF_B = EGMF_B  
        self.EGMF_LB = EGMF_LB
        if simu_path == "":
            self.path = make_path_dir(redshift,Emin,Emax,EBLmodel,EGMF_B,EGMF_LB)
        else:
            self.path = simu_path
            
        if Emin == Emax: 
            self.Erange = [Emin_simu,Emax]
            self.monoE = True
        else:
            self.Erange = [Emin,Emax]
            self.monoE = False
       
    def Run(self,force=False,z_0=False):
        '''
        Check if the simulation already exists or not.
        Run it if necessary and load data.
        This second step is mandatory before doing anything else!
        If you want to reinitialize any selection, jet or source spectrum applied, just 
        rerun this step.
        * force -> will rerun even if the simulation already exist
        * z_0 -> if the simulation has been designed to record photons at z=0 and not on 
            Earth. Require to recompile the program.  
        '''
        current_path = getcwd() 
        chdir(CECsi_PATH)
         
        if not path.exists(self.path) or force:
            print("Simulation does not exist", self.path)
            print("launch ", ctime())
            if not path.exists(self.path):
                makedirs(self.path) 
            write_preprocessing_f95(self.EBL)
            write_input_parameters_f95(self.redshift,self.Emin,self.Emax,self.Erange[0],self.EGMF_B,self.EGMF_LB,OMP_num_threads=5)

            subprocess.call("make")   
            subprocess.call("./CECsi.exe")   
            # copy results from simulation
            if path.isfile("temp/results.dat"):
                shutil.copy("temp/results.dat",self.path)
            if path.isfile("temp/results_z_0.dat"):
                shutil.copy("temp/results_z_0.dat",self.path)
            shutil.copy("temp/profile.dat",self.path)
            shutil.copy("temp/Compton_Errors.dat",self.path)
            shutil.copy("temp/input_parameters.f95",self.path)
            shutil.copy("temp/preprocessing.f95",self.path)
            if path.isfile("temp/leptons_loosed.dat"):
                shutil.copy("temp/leptons_loosed.dat",self.path)
            print("\n")

        self.LoadData(z_0)

        chdir(current_path)

    def LoadData(self,z_0=False):
        # Data contents --------------------------------------------------
        # 0 - generation
        # 1 - weight
        # 2 - energy [GeV]
        # 3 - time delay [s]
        # 4,5,6 - theta_pos, theta_dir, phi_dir [rad]
        # 7 - source photon energy [GeV]
        # 8 - charge (for specific use)
        #-----------------------------------------------------------------
        if z_0:
            self.data = loadtxt(self.path+resultsFile_z_0,unpack=True)
        else:
            self.data = loadtxt(self.path+resultsFile,unpack=True)

        if self.monoE:
            self.data[1] /= 10000 
        else:
            self.data[1] /= 50000 
 

    def Apply_source_spectrum(self,Es,dNdEs):
        '''
        Allow to apply a source spectrum to simulation data.
        Es and dNdEs are tables giving source energy and dN/dEs corresponding.
        Make sure that Es table is larger or equal to the Emin, Emax of the source.
        If the source is mono-energetic (Emin = Emax), this won't apply.
        '''
        if (min(Es) <= self.Erange[0]) and (max(Es) >= self.Erange[1]): 
            self.data[0] *= 10**(interp(log10(self.data[7]),log10(Es),log10(dNdEs)))
        else:
            print("Given Es range is too short compare to the simulated source energy range: [%e,%e]"%(self.Erange[0],self.Erange[1]))

    def Apply_jet(self,tjet=180.,tobs=0.,image=False):
        '''
        Allow to apply a jet to simulation data.
        * tjet [deg.] -> jet opening angle
        * tobs [deg.] -> jet misalignement with respect to the line of sight
        image [bool. optional] -> must be set to True if you to draw picture otherwise it 
        slow down the treatment of data
        '''
        if image:
            Nphi=int(150*(tobs/45.+1)**2)
            saveRes = []  
            phi_dir = self.data[6].copy()
            for phi_dj in (2*arange(float(Nphi))/Nphi-1)*pi:
                # replace phi_dir in source frame by phi_dir in observer frame
                self.data[6] = phi_dj 
                cos_theta_e = cos(tobs*degre)*cos(self.data[4])+sin(tobs*degre)*sin(self.data[4])*cos(phi_dj-phi_dir)
                cond = (arccos(cos_theta_e)/degre < tjet) # jet selection
                if saveRes == []:
                    saveRes = self.data[:,cond]
                else:
                    saveRes = append(saveRes,self.data[:,cond],axis=1)
                    self.data = saveRes    
                    self.data[1] /= Nphi 

        else:
            if tjet != 180:
                # in the case of a conic jet, some events can be rejected without computed 
                # all the phi_d cases; 
                #if cos(theta_obs - theta_pos^code) < cos(tjet) -> reject
                cond = cos(tobs*degre-self.data[4,:]) >= cos(tjet*degre)
                self.data = self.data[:,cond]
                cond = (cos(tobs*degre+self.data[4,:]) < cos(tjet*degre)) & (cos(tobs*degre-self.data[4,:]) > cos(tjet*degre))
                self.data[1,cond] *= arccos((cos(tjet*degre)-cos(tobs*degre)*cos(self.data[4,cond]))/(sin(tobs*degre)*sin(self.data[4,cond])))/pi

        # apply opening angle 
        cond = self.data[5] < fov*degre
        self.data = self.data[:,cond]

    def Data_selection(self,Erange=[1e-3,1e6],delayrange=[-1e8,1e30],fov=180.,generation='all',charge=0):
        '''
        Apply a simple selection of the simulation.
        Erange [GeV] -> energy min and max 
        dealyrange [s] -> time range 
        fov [deg.] -> field of view or psf
        generation -> photon generation desired
        charge -> 0 = photon, -1 = electron, 1 = positron
        '''
        if charge == 0:
            try:
                gen_select = int(generation)
                cond = (self.data[0] == gen_select*2) # only photons 
            except:
                # select photons only
                cond = (self.data[0]%2==0) 
        else:
            cond = (self.data[8]==charge)  
            try:
                gen_select = int(generation)
                cond = cond & (self.data[0] == gen_select*2+1) # only leptons
            except:
                # select leptons only
                cond = cond & (self.data[0]%2!=0) 
        # energy selection, time selection and angle selection
        cond = cond & (self.data[2] >= Erange[0])     & (self.data[2] <= Erange[1])  
        cond = cond & (self.data[3] >= delayrange[0]) & (self.data[3] <= delayrange[1])
        cond = cond & (self.data[4] <= fov*degre) 
        self.data = self.data[:,cond]

    def EventsList(self):
        '''
        Return events dataset. 
        Can be used at any time after running the simulation (CECsimulation.Run())
        To draw images, use CECsimulation.Apply_jet with option "image=True" before.

        Output are:
        [0] weight
        [1] energy [GeV]
        [2] time delay [s]
        [3] arrival angle [degre] 
        [4,5] theta, phi (arrival polar and azimuthal coordinates) [degre]
        [6] generation 
        '''
        return self.data[1], self.data[2], self.data[3], self.data[5]*rad, self.data[6]*rad, self.data[7]*rad, self.data[0]/2

    def Distribution(self,datatype,nbBins=100,x_range=[]):
        '''
        compute the 
        [datatype = "energy"] energy spectrum 
        or [datatype = "angle"] the angular distribution 
        or [datatype = "time"] the time distribution.

        NbBins -> number of bins
        x_range -> data range to compute the distribution. By default it is the min and
        max of the dataset. 
        '''
        if datatype == "energy":
            x = self.data[2]
        if datatype == "angle":
            x = self.data[5]*rad
        if datatype == "time":
            x = self.data[3]/yr
        weight = self.data[1]

        if x_range==[]:
            x_range=[min(x[x>0]),max(x)]
        dN,xi = histogram(log10(x[x>0]),nbBins,range=log10(x_range),weights=weight[x>0])
        xi = 10**xi
        xc = (xi[1:nbBins+1]+xi[0:nbBins])/2
        dx = xi[1:nbBins+1]-xi[0:nbBins]
        return xc, dN/dx

    def does_it_exist(self):
        '''
        Test if the simulation has already been runned.
        '''
        return path.exists(self.path)
