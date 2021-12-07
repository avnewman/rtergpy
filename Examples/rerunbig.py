#!/usr/bin/env python3

"""
Uses the std output from GT's searchCMT algo to pull and process multiple events.
Andrew Newman Fri Oct 29 13:17:32 EDT 2021

It took 10.25 minutes to process 17 events including pulling data from NEIC (these were done in series)
"""
from rtergpy.run import defaults, event, etime2name, src2ergs
from obspy import UTCDateTime
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('Agg')   # needed to plot without using the X-session

Defaults=defaults()
Event=event()
#Defaults.waveparams=[ [[1/300,2],[0.5,2]], [-60,300] ,1 ]  # change frequencies
Defaults.waveparams=[ [[1/300,0.5],[0.5,1]], [-60,600], 1 ]  # change frequencies and duration
Defaults.src='IRIS'
Event.newData=True # True is default
Event.ecount='00'
Event.iter="01"
edateold=""
MMIN=8.0
# events older than available in NEIC
CMTS=pd.read_csv('CMTS_IRIS.txt', sep='\s+', comment="#")  # any amount of whitespace
#Defaults.src='NEIC'
#CMTS=pd.read_csv('CMTS.txt', sep='\s+', comment="#")  # any amount of whitespace
for index, EQ in CMTS.iterrows():
    eloc = [EQ.LAT,EQ.LONG,EQ.DEPTH] 
    #print(EQ.DATE,EQ.TIME)
    year,mo,dd = EQ.DATE.split('/')
    hh,mn,sec = EQ.TIME.split(':')
    etime=(UTCDateTime(int(year),int(mo),int(dd),int(hh),int(mn),float(sec)))
    # iterate ecount
    if EQ.DATE == edateold:
        Event.ecount=str(int(Event.ecount)+1).zfill(2)
    else:
        Event.ecount='00'
    edateold=EQ.DATE
    Event.eventname=etime2name(etime,ecount=Event.ecount)
    Event.origin=[eloc,etime]
    Event.focmech=[EQ.STK, EQ.DP, EQ.RKE] # phi,delta,lmbda

    if EQ.Mw >= MMIN:    # only run on really large events
    	print("\n\n"+Event.eventname+" Mw= "+str(EQ.Mw)+" ===============================")
    	#try:
    	src2ergs(Defaults=Defaults,Event=Event)  # need to export run output in a coherent way
    	plt.close('all')  # they don't close themselves
    #except:
    #    print("ERROR: running on "+Event.eventname+" failed!!!!\n\n")
