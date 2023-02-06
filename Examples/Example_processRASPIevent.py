#!/usr/bin/env python3

"""
program to get event-based waveform data around theoretical P-arrival time
from the NEIC server, and including instrument response removal from IRIS
to the traces, the following informaiton is added:
Stations currently default to the GSN network
 Andrew Newman  Fri Oct 29 12:35:11 EDT 2021
"""
from rtergpy.run import defaults, event, etime2name, src2ergs
from obspy import UTCDateTime

Defaults=defaults()
Event=event()

Defaults.src="RASPISHAKE"
Defaults.network="AM"
Defaults.chan="EHZ"
# 2023 Turkey EQ (Feb 6, 2023 - mainshock)
eloc = [37.17,37.03,17.9]
etime= UTCDateTime(2023,2,6,1,17,35) 
Event.ecount='00'
#Event.newData=False   # use already downloaded data
Event.newData=True
Event.eventname=etime2name(etime,ecount=Event.ecount)
Event.origin=[eloc,etime]
Event.focmech=[318, 89, -179] # phi,delta,lmbda

print(Event.eventname)
src2ergs(Defaults=Defaults,Event=Event)  # need to export run output in a coherent way
