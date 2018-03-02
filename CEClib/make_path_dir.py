def make_path_dir(redshift,Emin,Emax,EBL,EGMF_B,EGMF_LB):
   '''
   Determine the subfolder path where the simulation is or will be.
   '''
   simu_dir = "MCsimulations_DB/%s/z=%1.9f/"%(EBL,redshift) 
   LB = str(int(EGMF_LB)) 
   if EGMF_LB <= 1e-6:
      LB = "%1.6f"%EGMF_LB
   elif EGMF_LB <= 1e-5:
      LB = "%1.5f"%EGMF_LB
   elif EGMF_LB < 1:
      LB = str(EGMF_LB) 
   if Emin != Emax:
      simu_dir += "EGMF=%1.0eG_"%(EGMF_B)+LB+"Mpc_Emax=%1.0fTeV"%(Emax)
   else:
      simu_dir += "EGMF=%1.0eG_"%(EGMF_B)+LB+"Mpc_Emax=%1.0fTeV_mono"%(Emax)
   return simu_dir      
