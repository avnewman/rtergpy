"""
Algorithm to use Obspy's metadata to pull response and other metadata.  Returns an apporpiate
inventory class.
A.V. Newman Mon Jul 26 15:26:35 EDT 2021
"""
from obspy.clients.fdsn import Client as fdsnClient
from obspy import UTCDateTime

def get_respinv(network,eloc,etime,rads,chan):
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
