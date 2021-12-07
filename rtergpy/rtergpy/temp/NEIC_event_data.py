#!/usr/bin/env python3

# program to get event-based waveform data around theoretical P-arrival time
# from the NEIC server, and including instrument response removal from IRIS
# to the traces, the following informaiton is added:
#      stats.phasePtime  (using IASP91)
#      stats.coordinates ({'latitude': slat, 'longitude': slon})
#      stats.distance    (in meters around great circle)
# Stations currently default to the GSN network
#  Andrew Newman  Fri Jul 23 18:56:10 UTC 2021

from obspy.clients.fdsn import Client as fdsnClient
from obspy import UTCDateTime
from obspy.taup import TauPyModel
from obspy.core.stream import Stream
model = TauPyModel(model="iasp91")
from obspy.geodetics.base import locations2degrees as l2d
from obspy.geodetics.base import degrees2kilometers as d2km
from obspy.clients.neic import Client as nClient
nclient=nClient()
from compress_pickle import dump,load  # reading/writing compressed pickles

### params
### need to determine best method for program arch.  May best be done as all function calls
# event
eloc = [7.39,-82.5,33]
etime= UTCDateTime(2021,7,21,21,15,15)   # recent mid-6 in Panama
eloc = [13.80,120.50,104]
etime= UTCDateTime(2021,7,23,20,49,00)   # recent mid-6 in Panama

## general
network = "CU,GT,IC,II,IU" # GSN network
rads=[25,80]
chan="BHZ"
prePtime=60
postPtime=300

def get_inventory(network,eloc,etime,rads,chan):
    fclient = fdsnClient()
    now=UTCDateTime()

    elat,elon,edep = eloc
    minrad,maxrad = rads 
    
    inventory = fclient.get_stations(network = network, 
            latitude = elat, longitude = elon, 
            minradius = minrad, maxradius = maxrad, 
            starttime=etime-86400, endtime=etime,
            channel=chan,
            #location="00",
            matchtimeseries=True,
            #filename="test.txt", format="text",  # cannot be used with response
            #filename="test.xml", format="xml",
            level="response"
            )
    return inventory  # this is an Obspy type

def theorPtime(eloc,etime,sloc):
    elat,elon,edep=eloc
    slat,slon,sheight=sloc
    arrivals=model.get_travel_times_geo(
        source_depth_in_km=edep, source_latitude_in_deg=elat, source_longitude_in_deg=elon,
        receiver_latitude_in_deg=slat, receiver_longitude_in_deg=slon,
        phase_list="P")
    return etime + arrivals[0].time

def loc2km(eloc,sloc):
    elat,elon,edep=eloc
    slat,slon,sheight=sloc
    distdeg=d2km(l2d(elat,elon,slat,slon))
    return distdeg
    
# time run
t1 = UTCDateTime()

inventory = get_inventory(network,eloc,etime,rads,chan)
#print(inventory)
st=Stream()
st.clear()  # clear just in case

# run on all channels in inventory
for chan in inventory.get_contents().get("channels"):
    slat,slon,sz,sldep=inventory.get_coordinates(chan).values()
    sheight=sz-sldep
    sloc=slat,slon,sheight
    Ptime=theorPtime(eloc,etime,sloc)
    neti,stati,loci,chani=chan.split(".")
    stlocal=nclient.get_waveforms(neti,stati,loci,chani,Ptime-prePtime,Ptime+postPtime)
    if stlocal:
        # add station coordinates
        try:
            stlocal[0].stats.coordinates= {'latitude': slat, 'longitude': slon}
        except:
            print ("Channel ", chan, " did not add event locations")
        # add distance
        try:
            stlocal[0].stats.distance=loc2km(eloc,sloc)*1000;  # distance should be reported in meters!
        except:
            print ("Channel ", chan, " did not calc distance")
        # add theoretical P-time
        try:
            stlocal[0].stats.phasePtime=Ptime;  # UTC time of P-arrival
        except:
            print ("Channel ", chan, " Ptime not added")
        st+=stlocal
        print("added: ",chan)
    else:
        print("- skip - :", chan)
        
st_raw=st.copy()  # create backup

# process data
st.attach_response(inventory) # include response information
st.detrend(type='demean')
st.taper(0.05)
st.remove_response()
st.filter("bandpass", freqmin=0.01, freqmax=2)

print("\nRuntime: ", UTCDateTime()-t1, " sec.")


# save results
#st.write("processed.mseed", format="MSEED")  # these may not save some components (e.g. phasePtime)
#st_raw.write("raw.mseed", format="MSEED")
st.write("processed.pkl", format="PICKLE")
st_raw.write("raw.pkl", format="PICKLE")

# compression is fast and about 2x smalelr
dump(st,"processed.pkl.gz", compression="gzip")
dump(st_raw,"raw.pkl.gz", compression="gzip")

# test 
exit(0)
