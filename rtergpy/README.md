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
pip install update .
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