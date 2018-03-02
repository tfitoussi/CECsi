def write_preprocessing_f95(EBL):
   '''
   Write preprocessing.f95 file with the correct EBl model
   '''
   F = open('temp/preprocessing.f95','w')
   F.write("! Choose EBL model to use (ONLY ONE MUST UNCOMMENT!) =============\n")
   F.write("!#define EBL_model_Dominguez\n")
   F.write("!#define EBL_model_Finke\n")
   F.write("!#define EBL_model_Franceschini\n")
   F.write("!#define EBL_model_Gilmore\n")
   F.write("!#define EBL_model_Best_fit\n")
   F.write("!#define EBL_model_Lower_limit\n")
   F.write("#define EBL_model_%s\n"%EBL)
   F.write("\n")
   F.write("! Select EGMF model (ONLY ONE MUST UNCOMMENT!) ====================\n")
   F.write("!#define B_uniform\n")
   F.write("#define B_pseudo_random\n")
   F.write("!#define B_random\n")
   F.write("!#define B_turbulent\n")
   F.write("\n")
   F.write("! Approximations =================================================\n")
   F.write("! here activate or deactivate separately each approximation in the code\n")
   F.write("! To deactivate just add ! before #define ... \n")
   F.write("!#define approx_motionless_leptons\n")
   F.write("!#define approx_Ee\n")
   F.write("!#define approx_Eic\n")
   F.write("!#define approx_lambda_gg\n")
   F.write("!#define approx_lambda_ic\n")
   F.write("\n")
   F.write("! here add an output file which contains: leptons deflection angle,\n")
   F.write("! energy and d_ic\n")
   F.write("!#define file_lepton_deflection\n")
   F.write("!#define file_positrons\n")
   F.write("!#define file_cascade_traj\n")
   F.write("\n")
   F.write("! to set the limit to z=0 and not when the particle reachs the Earth\n")
   F.write("!#define z_0_limit\n")
   F.close()

def write_input_parameters_f95(redshift,Emin,Emax,Emin_sel,EGMF_B,EGMF_LB,OMP_num_threads=5):
   '''
   Write input_parameters.f95 file with the correct parameters (redshift, energy range, EGMF)
   N.B.: number of photons launched has been set arbitrarily to 
      10000 if the source is mono-energetic (Emin = Emax)
      50000 else
   it should give a sufficient statistic. 
   '''
   F = open('temp/input_parameters.f95','w')
   F.write("!==== Source ====\n")
   F.write("&source\n")
   F.write("   emission_redshift = %f,\n"%redshift) 
   F.write("   Emin_source = %f, !TeV\n"%Emin)
   F.write("   Emax_source = %f, !TeV\n"%Emax)        
   F.write("   Nature_source = 0,\n")            
   if Emin != Emax:
      Nmax = 50000
   else:
      Nmax = 10000
   F.write("   nmax = %d,\n"%Nmax)                
   F.write("&end\n")
   F.write("\n")
   F.write("!==== EGMF ====\n")
   F.write("&EGMF_param\n")
   F.write("   intEGMF = %1.0e, !G\n"%EGMF_B)
   F.write("   lambdaEGMF = %1.0e, !Mpc\n"%EGMF_LB)
   F.write("   Nm = 50,\n")     
   F.write("&end\n")
   F.write("\n")
   F.write("!==== Simulation parameters ====\n")
   F.write("&simulation\n")
   F.write("   OMP_num_threads = %d,\n"%OMP_num_threads) 
   F.write("   Ethreshold = %1.0e, !GeV\n"%Emin_sel)        
   F.write("   Comptonthreshold =0.005,\n") 
   F.write("   alphaS = 0.6,\n")
   F.write("   seed = 1235789,\n")
   F.write("&end\n")
   F.close()
