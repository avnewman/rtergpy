#!/bin/env python3
from rtergpy.run import mergeResults
import pandas as pd

runs=mergeResults(iteration='01')
CMTSIRIS=pd.read_csv('CMTS_IRIS.txt', sep='\s+', comment="#")
CMTSNEIC=pd.read_csv('CMTS_NEIC.txt', sep='\s+', comment="#")
CMTS=pd.concat([CMTSIRIS,CMTSNEIC], ignore_index=True)
runs=runs.assign(Mom=CMTS['Mom'])
runs=runs.assign(Mw=CMTS['Mw'])
runs=runs.assign(CMTNAME=CMTS['CMTNAME'])
runs.to_pickle("Results_2000-late2021.01.pkl")
print(runs)
