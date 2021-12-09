#!/usr/bin/env python3
from rtergpy.run import mergeResults
import pandas as pd

# will look in default location
runs=mergeResults(iteration='00')
CMTS=pd.read_csv('CMTS_NEIC.txt', sep='\s+', comment="#")
runs=runs.assign(Mom=CMTS['Mom'])
runs=runs.assign(Mw=CMTS['Mw'])
runs=runs.assign(CMTNAME=CMTS['CMTNAME'])
runs.to_pickle("Combined_Results.pkl")
print(runs)
