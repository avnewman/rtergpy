"""
Algorithm to use taup event, and station location info to calculate theoretical P-time
A.V. Newman Mon Jul 26 15:26:35 EDT 2021
"""
from obspy.taup import TauPyModel
model = TauPyModel(model="iasp91")

def theorPtime(eloc,etime,sloc):
    elat,elon,edep=eloc
    slat,slon,sheight=sloc
    arrivals=model.get_travel_times_geo(
        source_depth_in_km=edep, source_latitude_in_deg=elat, source_longitude_in_deg=elon,
        receiver_latitude_in_deg=slat, receiver_longitude_in_deg=slon,
        phase_list="P")
    return etime + arrivals[0].time
