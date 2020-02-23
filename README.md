# Cosmological Electromagnetic Cascades simulation (CECsi)


This program is a Monte-Carlo simulation of electromagnetic Cosmological cascades. A cosmological cascade is a simple process. High energy photons interact with low energy photons (mostly from the extragalactic background) to produce a pair of electron and positron. Each of the leptons has an energy equivalent to half of the initial photon energy. Then the lepton could interact via inverse Compton scattering with low energy photons (mostly form the Cosmological Microwave Background). This produce a photon with a higher energy (GeV) and the lepton loose the equivalent energy. Though each particle could interact again, respectively through pair production and inverse Compton scattering, or reach the detector.

** This code is an updated version of the one presented in this [article](http://adsabs.harvard.edu/cgi-bin/basic_connect?qsearch=fitoussi+2017&version=1) because it appears that the simulation code is not easy to use by itself and a lot of options have been add with time. This new version aims to be easier too used. The old code can still be found in the [old directory](https://gitlab.com/tfitoussi/cascade-simulation.git) but will most probably not be updated anymore. **

## Requirements

This program was written with the objective to be easy to use, as fast as possible and describing in the most realistic way the process without analytic simplification. The MC code is written in `Fortran 95`. It has been successfully compiled with `gfortran` under Linux, Mac and Solaris using `gmake`. Moreover it is parallelized with `OpenMP` which is embedded most of the fortran compilers.

The simulation returns series of events which need to be post-treated to generate spectra, images or timing delays. Post-treatment scripts are written in `python 3.5` and requires `Numpy/Scipy` + `Matplotlib` to draw figures. 

These programs are distributed under the General Public License V3. Feel free to use it, improve it and redistribute it.

**Disclaimer: this program comes without any guarantees. Beware of errors and use common sense and critical mind interpreting results!** 

Summary:

* Fortran compiler (gfortran is prefered with version higher than or equal to 4.8.5)   

* `make` (gmake is prefered)

* `Python 3.5`

* [`Numpy/Scipy` + `Matplotlib`](http://www.scipy.org/install.html) 

* `ipython` and `astropy` are not mandatory but recommanded


## Installation 

We strongly recommand to install CECsi in a python virtualenv. It is easier to manage and install python libraries. Moreover if you are working on a machine on which you don't have the administration right this is the easiest way.

If `virtualenv` is not install it or download it:

```bash
wget https://github.com/pypa/virtualenv/archive/develop.zip
unzip develop.zip
python virtualenv-develop/virtualenv.py $CRPROPA_DIR
```

### Setup the virtualenv

If you already install (e.g.) CRPropa in a virtualenv, you can install CECsi in it to work alongside with CRPropa. But you can still install CECsi in a separate virtualenv as far as CECsi is an independant software.

1. Create the virtualenv
```bash
export VENVDIR=$HOME"/.virtualenvs/venv"
mkdir -p $VENVDIR
virtualenv $VENVDIR
```
2. Activate the virtualenv
```bash
source $VENVDIR"/bin/activate"
```
3. You may want to add the alias to your `.bashrc`: add this two lines in your `.bashrc` file to launch the virtualenv with the command `cecsi`
```bash
export VENVDIR=$HOME"/.virtualenvs/venv"
alias venv="source $VENVDIR'/bin/activate'"
```

### Download CECsi (with or without virtualenv)

4. Download the code
```bash
cd $VENVDIR
git clone https://gitlab.com/tfitoussi/CECsi.git
```

### Finish to setup the virtualenv

5. Add CECsi path to the virtualenv python path 
```bash
cd CECsi
pwd > ../lib/python2.7/site-packages/CECsi.pth 
```

6. Install python dependencies
```bash
make pip
```

### Before launching any simulations (with or without virtualenv)

7. Create needed directories
```bash
make initiate
```

Congratulations! CECsi is now ready to use. 


## Using CECsi code

Some example on how-to use CECsi has been put in the Examples directory. Some explanations are writen inside. If you need help, do not hesitate to contact me. 
