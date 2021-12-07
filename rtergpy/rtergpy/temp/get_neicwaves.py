"""
Algorithm to use Obspy's metadata from get_respinv (fsdn) to pull data from
the NEIC server around the P-wave theoretical time.  
Returns an streams with added station location, event distance metadata.

A.V. Newman Mon Jul 26 15:26:35 EDT 2021
"""

from obspy.core.stream import Stream
from obspy.clients.neic import Client as nClient
from rtergpy import theorPtime, loc2km
nclient=nClient()

def get_neicwaves(inventory,eloc,etime,preP,postP):
    st=Stream()

    # run on all channels in inventory
    for chan in inventory.get_contents().get("channels"):
        slat,slon,sz,sldep=inventory.get_coordinates(chan).values()
        sheight=sz-sldep
        sloc=slat,slon,sheight
        Ptime=theorPtime(eloc,etime,sloc)
        neti,stati,loci,chani=chan.split(".")
        stlocal=nclient.get_waveforms(neti,stati,loci,chani,Ptime-preP,Ptime+postP)
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
        
    st.attach_response(inventory) # include response information
    return st 