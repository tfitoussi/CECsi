# Cosmological Electromagnetic Cascades simulation (CECsi)


This program is a Monte-Carlo simulation of electromagnetic Cosmological cascades. A cosmological cascade is a simple process. High energy photons interact with low energy photons (mostly from the extragalactic background) to produce a pair of electron and positron. Each of the leptons has an energy equivalent to half of the initial photon energy. Then the lepton could interact via inverse Compton scattering with low energy photons (mostly form the Cosmological Microwave Background). This produce a photon with a higher energy (GeV) and the lepton loose the equivalent energy. Though each particle could interact again, respectively through pair production and inverse Compton scattering, or reach the detector.

**Because it appears that the simulation code is not easy to use by itself and requires to use post-treatment scripts which were stored in a different [git repository](https://gitlab.com/tfitoussi/simulation-analysis.git) has say in the [article](http://adsabs.harvard.edu/cgi-bin/basic_connect?qsearch=fitoussi+2017&version=1). In a way of simplification, everything has been merge here and an additional script has been add to help running thing out of the box.**

## Requirements

This program was written with the objective to be easy to use, as fast as possible and describing in the most realistic way the process without analytic simplification. The MC code is written in `Fortran 95`. It has been successfully compiled with `gfortran` under Linux, Mac and Solaris using `gmake`. Moreover it is parallelized with `OpenMP` which is embedded most of the fortran compilers.

The simulation returns series of events which need to be post-treated to generate spectra, images or timing delays. Post-treatment scripts are written in `python 3.5` and requires `Numpy/Scipy` + `Matplotlib` to draw figures. 

These programs are distributed under the General Public License V3. Feel free to use it, improve it and redistribute it.

**Disclaimer: this program comes without any guarantees. Beware of errors and use common sense and critical mind interpreting results!** 

Summary:

* Fortran compiler (gfortran is prefered with version higher than or equal to 4.8.5)   

* `make` (gmake is prefered)

* `Python 3.5` (`Python 2.7` should work too)

* [`Numpy/Scipy` + `Matplotlib`](http://www.scipy.org/install.html) 


## First use

### Installation in python Virtualenv (Recommended)

1. Create the virtualenv

   export CECSI_DIR=$HOME"/.virtualenvs/CECsi"

   mkdir -p $CECSI_DIR

   virtualenv $CECSI_DIR

2. Activate the virtualenv

   source $CECSI_DIR"/bin/activate"

3. add the alias to your `.bashrc`: add this two lines in your `.bashrc` file to launch the virtualenv with the command `cecsi`

   export CECSI_DIR=$HOME"/.virtualenvs/CECsi"
   alias cecsi="source $CECSI_DIR'/bin/activate'"

4. Download and install the code

   cd $CECSI_DIR

   git clone git@gitlab.com:tfitoussi/CECsi.git

5. Add CECsi path to the virtualenv python path

   cd CECsi

	pwd > ../lib/python2.7/site-packages/CECsi.pth 

### Before launching any simulations (with or without virtualenv)

6. `Create needed directories

   mkdir -p "temp" "Modules" "MCsimulations_DB"

7.  Check that the MC simulation code is successfully compiled

   make


## Using CECsi code

Some example on how-to use CECsi has been put in the Examples directory. All explanations are writen inside.


