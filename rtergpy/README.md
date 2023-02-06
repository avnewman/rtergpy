# About

This github provides `rtergpy`; a set of python codes for the calculation of the radiated seismic energy as determined at the approximated termination of earthquake rupture.   

The code is modified signficantly from the original FORTRAN development for the Energy/Moment Ratio parameter (Theta parameter) as described in [Newman and Okal (1998)](http://geophysics.eas.gatech.edu/people/anewman/research/papers/Newman_Okal_JGR_1998.pdf).  
The code now calculates the per-second energy change as determined per station ([Convers and Newman, (2011)](http://geophysics.eas.gatech.edu/people/anewman/research/papers/Convers_Newman_JGR_2011.pdf)), and the total energy
radiated at the maximum in the Time-Averaged Cumulative Energy Rate (TACER) as decribed in [Convers and Newman (2013)](http://geophysics.eas.gatech.edu/people/anewman/research/papers/Convers_Newman_GRL_2013.pdf).


A paper describing this new code is forthcoming.

## Install 
### Create a usable Conda environment for working with the package
```
# create with a bunch of packages
conda create --name rtergpy python=3.11 ipython ipykernel 
```
### Once you're in the proper environment, for me it is done by:
``` 
conda activate rtergpy
```
### First install gmt using conda
```
conda install gmt
```
### Install remaining external packages via pip
```
pip install pygmt compress_pickle matplotlib tqdm obspy pandas
```

### Install the `rtergpy` pagkage by going into this directory `rtergpy` (directory with `setup.py`) and using pip:

#### Install for general user:
```
pip install .
# alternatively, install from another directory, but point to that directory
pip install '/path/to/setup/.py/file'
```
If there are any additonal dependencies needed, you may need to install these first.

To update your install, after downloading (e.g. through git pull)
```
# remember to be in the `setup.py` directory
pip install update .
```
#### Install for developer:
Alternatively, if you are planning on modifying the code for your own development, it would likely be better to run the below, to allow updates to code to be automatically loaded

Run in 'Edit' Mode
```
# Allowing edits to programs to automatically become loaded
pip install -e .
```
----
## Setup
Before running any code, copy `rtergpy.conf` into your base directory and make any changes to that version for
code and event directory locations, as well as any other parameters that you want to change for your default runs. 

 `cp rtergpy.conf ~/` 
 
 > Note that program will check `/etc/`, then `~/`, finally the local run directory (`.`) for config files, overriding configuration details along the way. Thus, directory specific configurations will always supercede prior configs.
 
----
## Run
Look in the `Examples` and `notebooks` for how to process. 
