#!/usr/bin/env python3
"""
Example script to calculate radiated energy from seismo-geodetic data using RTERGPY,
adapted for locally stored velocity data from GNSS.
"""

import os
import pandas as pd
import numpy as np
from obspy import UTCDateTime
from pathlib import Path

from rtergpy.run import defaults, event, etime2name, src2ergs

# === Set your event name ===
event_name = "Tohoku_20110311" 

# === Define base paths ===
base_path = "/Users/hkunwer/Documents/research/EQenergy/seismogeodetic/data"
event_path = os.path.join(base_path, event_name)
event_metadata = os.path.join(event_path, "metadata")
event_mseed = os.path.join(event_path, "kvel", "mseed")
meta_filename = f"{event_name}.kvel.meta.csv"
df_path = os.path.join(event_metadata, meta_filename)

# === Load event metadata from EarthquakesTable ===
earthquake_data = pd.read_csv("EarthquakesTable.csv")
event_data = earthquake_data[earthquake_data["Event_name"] == event_name].iloc[0]

origin_time = UTCDateTime(event_data["Origin_time(UTC)"])
year = origin_time.year
sensor_type = "kvel"

# === Create output directory ===
output_dir = Path("events") / str(year) / event_name / sensor_type
output_dir.mkdir(parents=True, exist_ok=True)
print(f"Output will be saved in: {output_dir}")

# === Initialize RTERGPY objects ===
Defaults = defaults()
Event = event()

Defaults.src = "Seismo-geodetic"
Defaults.manual = True
Defaults.remove_response = False
Defaults.savedir = output_dir

# === Define event properties ===
eloc = [event_data['Latitude(N)'], event_data['Longitude(E)'], abs(event_data['Depth(km)'])]
etime = UTCDateTime(event_data['Origin_time(UTC)'])

Event.eloc = eloc
Event.etime = etime
Event.ecount = '00'
Event.eventname = etime2name(etime, ecount=Event.ecount)
Event.origin = [eloc, etime]
Event.iter = "GEOD"
Event.eventdir = event_path
Event.is_seismogeodetic = True  # Use seismoGNSS path
Event.mseed = event_mseed
Event.dfpath = df_path

# === Extract and assign focal mechanism ===
try:
    foc_str = event_data["Focal_mechanism"]
    strike, dip, rake = map(float, foc_str.strip().split("/"))
    Event.focmech = [strike, dip, rake]
    print(f"✅ Focal mechanism set: {Event.focmech}")
except Exception as e:
    print(f"⚠️ Failed to parse focal mechanism for {event_name}: {e}")
    Event.focmech = None

# === Run radiated energy calculation ===
print(f"Running energy calculation for {Event.eventname}...")
src2ergs(Defaults=Defaults, Event=Event)
