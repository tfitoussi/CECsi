from numpy import size, chararray, array, where, log10
from read import readxml

def diff_simu(files):
   label = chararray(size(files),itemsize=100)
   simu0 = array(readxml("Results/"+files[0]+"/info.xml"))
   for i in range(1,size(files)):
      simu = array(readxml("Results/"+files[i]+"/info.xml"))
      Ni = where(simu!=simu0)[0]
      if Ni.size == 0:
         label[0] = simu0[0][0]
         label[i] = simu[0][0]
      else:
         label[0]=""
         label[i]=""
              
         if 1 in Ni: # redshift 
            label[0] += "z=%1.2f "%float(simu0[1][0])
            label[i] += "z=%1.2f "%float(simu[1][0])
            if Ni.size >2:
               label[0] += ", "
               label[i] += ", "
         if 3 in Ni: # Emax
            label[0] += "$E_{max}$=%1.0f TeV "%float(simu0[3][0])
            label[i] += "$E_{max}$=%1.0f TeV "%float(simu[3][0])
            if Ni.size >2:
               label[0] += ", "
               label[i] += ", "
         if (10 in Ni or 11 in Ni) and (3 not in Ni): # Emin_sel or Emax_sel
            label[0] += "[%1.1eGeV;%1.0eGeV] "%(float(simu0[10][0]),float(simu0[11][0]))
            label[i] += "[%1.1eGeV;%1.0eGeV] "%(float(simu[10][0]),float(simu[11][0]))
            if Ni.size >2:
               label[0] += ", "
               label[i] += ", "
         if 4 in Ni: # spectrum
            label[0] += parse_spectrum(simu0[4][0])
            label[i] += parse_spectrum(simu[4][0])
            if Ni.size >2:
               label[0] += ", "
               label[i] += ", "
         if 5 in Ni or 6 in Ni: # tjet or tobs
            label[0] += "$\\theta_{jet}$=%1.1f$^\circ$/$\\theta_{obs}$=%1.0f$^\circ$ " %(float(simu0[5][0]),float(simu0[6][0]))
            label[i] += "$\\theta_{jet}$=%1.1f$^\circ$/$\\theta_{obs}$=%1.0f$^\circ$ " %(float(simu[5][0]),float(simu[6][0]))
            if Ni.size >2:
               label[0] += ", "
               label[i] += ", "
         if 7 in Ni: # EBL
            label[0] += which_EBL(simu0[7][0])
            label[i] += which_EBL(simu[7][0])
            if Ni.size >2:
               label[0] += ", "
               label[i] += ", "
         if 8 in Ni or 9 in Ni: # EGMF B or LB  
            label[0] += "EGMF: 10$^{%1.0f}$G / %1.0eMpc "%(log10(float(simu0[8][0])),float(simu0[9][0]))
            label[i] += "EGMF: 10$^{%1.0f}$G / %1.0eMpc "%(log10(float(simu[8][0])),float(simu[9][0]))
            if Ni.size >2:
               label[0] += ", "
               label[i] += ", "
         if 12 in Ni: # fov
            label[0] += "fov=%1.1f$^\circ$ "%float(simu0[12][0])
            label[i] += "fov=%1.1f$^\circ$ "%float(simu[12][0])
   return label             
                

def which_EBL(EBL_model):
   if EBL_model == "Dominguez":
      return "Dominguez et Al (2011) "
   elif EBL_model  == "Finke":
      return "Finke et Al (2010) "
   elif EBL_model == "Franceschini":
      return "Franceschini et Al (2008) "
   elif EBL_model == "Gilmore":
      return "Gilmore et Al (2012) "
   elif EBL_model == "Best_fit":
      return "Kneiske et Al (2004) 'best fit' "
   elif EBL_model == "Lower_limit":
      return "Kneiske et Al (2004) 'lower limit' "

def parse_spectrum(spectrum):
   spectrum_attribute = spectrum.split(':')
   if spectrum_attribute[0]=="powerlaw":
      return "powerlaw ($\\Gamma$=%1.1f) "%float(spectrum_attribute[1]) 
   elif spectrum_attribute[0]=="logparabola":     
      return "log parabola ($E_{pivot}$=%1.1fGeV, $\\alpha$=%f, $\\beta$=%f) " %(float(spectrum_attribute[1]), float(spectrum_attribute[2]), float(spectrum_attribute[3]))
