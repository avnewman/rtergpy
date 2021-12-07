#!/usr/bin/env python3
"""
program uses currently processed stations and event data to calculate energy
at each station (or that's the plan)

Andrew Newman Wed Aug  4 16:12:23 EDT 2021
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from obspy import UTCDateTime
from compress_pickle import load  # reading/writing compressed pickles
from rtergpy.waveforms import wave2energy,gmeanCut,e2Me,loadwaves
import pandas as pd

cutoff=15
eloc = [13.80,120.50,104]
eventname='2021072300'  # data pulled py Example_rapidNEICeventdata.py
iter='00'
st,df=loadwaves(eventname,runiter=iter) # assumes iteration zero if not include
cutoff=15


#df=pd.DataFrame(dfline, dtype=object)
# run
#waveform window around P-wave
tstart=-5  # relative to P
tlength=70 # after tstart
tstep=1 # incriment
window=[tstart,tlength]
elat,elon,edepth=eloc
#additional event info
phi,delta,lmbda=[0, 20, 90]
foc=[phi,delta,lmbda]

#site conditions
pvel_site=7000.
rho_site=3000.

# run values
pmin=0.5; pmax=70;  # period range in seconds
fmin=1/pmax ; fmax=1/pmin; # freq range in Hz

t1 = UTCDateTime() # time run

E=np.empty((0,4))
for tr in st:
    thisE=wave2energy(tr,[[fmin,fmax],window],[eloc,foc],[pvel_site,rho_site])
    E=np.append(E,np.array([thisE]),axis=0)
 
#thisE=wave2energy(st,[[fmin,fmax],window],[eloc,foc],[pvel_site,rho_site])
now=UTCDateTime()
runtime=now-t1
print("\nRuntime: %.2f s\n" %runtime)

EgmeanTrue,ETrue=gmeanCut(E[:,1],cutoff=cutoff)
EgmeanEst,EEst=gmeanCut(E[:,0],cutoff=cutoff)
print("Estimated Energy =%.2e [Me=%.2f], Corrected Energy =%.2e [Me=%.2f]" %(EgmeanEst, e2Me(EgmeanEst),EgmeanTrue,e2Me(EgmeanTrue)) )
print("Using %d of %d observations" %(len(EEst),len(E)) )
plt.plot(ETrue,'bo',label="True Mech",alpha=0.3)
plt.plot(EEst,'r+',label="Est. Mech correction",alpha=0.3)
plt.yscale('log');
plt.legend();
plt.show();

#exit(0)
