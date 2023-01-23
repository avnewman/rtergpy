# About

This github provides`rtergpy`; a set of python codes for the calculation of the radiated seismic energy as determined at the approximated termination of earthquake rupture.   

The code is modified signficantly from the original FORTRAN development for the Energy/Moment Ratio parameter (Theta parameter) as described in [Newman and Okal (1998)](http://geophysics.eas.gatech.edu/people/anewman/research/papers/Newman_Okal_JGR_1998.pdf).  
The code now calculates the per-second energy change as determined per station ([Convers and Newman, (2011)](http://geophysics.eas.gatech.edu/people/anewman/research/papers/Convers_Newman_JGR_2011.pdf)), and the total energy
radiated at the maximum in the Time-Averaged Cumulative Energy Rate (TACER) as decribed in [Convers and Newman (2013)](http://geophysics.eas.gatech.edu/people/anewman/research/papers/Convers_Newman_GRL_2013.pdf).


A paper describing this new code is forthcoming.

## Install 
Once you're in the proper environment, for me it is done by:
``` 
conda activate rterg
```
Install by going into this directory `rtergpy` (directory with `setup.py`) and using pip:
```
pip install .
```
If there are any dependencies needed, you may need to install these first.

### Update
```
# remember to be in the `setup.py` directory
pip install update .
```

### Run in 'Edit' Mode
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
