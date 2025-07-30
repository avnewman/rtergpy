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

Defaults.src='IRIS'
# 2021 Haiti EQ *Mw 7.2)
# 21081402 2021/08/14 12:29:08  18.56  -73.55   10.0
#eloc = [18.56,-73.55,10] 
#etime= UTCDateTime(2021,8,14,12,29,8) 
eloc = [52.53, 160.165,20.7] 
etime= UTCDateTime(2025,7,29,23,24,50) 
Event.ecount='00'
Event.focmech=[198, 18,51] # phi,delta,lmbda
#Event.newData=False   # use already downloaded data
Event.newData=True
Event.eventname=etime2name(etime,ecount=Event.ecount)
Event.origin=[eloc,etime]
Event.Mi=8.8

print(Event.eventname)
src2ergs(Defaults=Defaults,Event=Event)  # need to export run output in a coherent way

