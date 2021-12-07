"""
Algorithm to use Obspy's geodetic distance calc for location in km
A.V. Newman Mon Jul 26 15:26:35 EDT 2021
"""
from obspy.geodetics.base import locations2degrees as l2d
from obspy.geodetics.base import degrees2kilometers as d2km

def loc2km(eloc,sloc):
    elat,elon,edep=eloc
    slat,slon,sheight=sloc
    distdeg=d2km(l2d(elat,elon,slat,slon))
    return distdeg