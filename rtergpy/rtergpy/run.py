from obspy import UTCDateTime
from tqdm import tqdm
import numpy as np
import pandas as pd
import os
import matplotlib
from configparser import ConfigParser

def etime2name(etime,ecount='00',**kwargs):
    """
    module converts UTCdateTime info (etime) and optionally count (ecount) into an eventname.
    If no ecount is supplied, code will assume ecount=0.
    """
    eventname=(str(etime.year)+str(etime.month).zfill(2)+str(etime.day).zfill(2)+str(ecount).zfill(2))
    return eventname

# Basic parameters for all processing

class defaults:
    def __init__(self):
        # read config file
        config = ConfigParser()
        # read etc, then base (~/), then local (overwriting along the way) 
        filelocations = ["/etc/rtergpy/", os.path.expanduser("~"), os.curdir ]
        try: 
            for dir in filelocations:
                file=str(os.path.join(dir,'rtergpy.conf'))
                config.read(file)
            DEFS = config["DEFAULTS"]
        except:
            print("ERROR: Could not find configuration file 'rtergpy.conf'")
            pass
       
        self.basedir=DEFS["basedir"]
        self.edirbase=DEFS["edirbase"]  
        self.libdir=DEFS["libdir"]

        self.network=DEFS["network"]
        self.stationrange=[ float(DEFS["statmin"]), float(DEFS["statmax"]) ] #[25,80]  # distance in degrees
        self.chan=DEFS["chan"]

        #site conditions
        pvel_site=float(DEFS["pvel_site"])
        rho_site=float(DEFS["rho_site"])
        self.siteparams=[pvel_site,rho_site]
        
        # frequency range and other waveform params
        f1min=float(DEFS["f1min"]) ; f1max=float(DEFS["f1max"]) # BB
        f2min=float(DEFS["f2min"]) ; f2max=float(DEFS["f2max"]) # HF 
        fbands=[[f1min,f1max],[f2min,f2max]]
        tstep=int(DEFS["tstep"]) # increment
        prePtime=int(DEFS["prePtime"])
        postPtime=int(DEFS["postPtime"])
        pwindow=[prePtime,postPtime]
        self.snr=float(DEFS["snr"])
        self.waveparams=[fbands,pwindow,tstep]
        self.resample=int(DEFS["resample"])  # samples per second
        self.smoothkern=int(DEFS["smoothkern"])  # kernel for gaussian smoothing of data for duration estimates (1/2 on each side)

        # earth params
        self.rearth=float(DEFS["rearth"]) # in meters
        self.qbc=float(DEFS["qbc"]) #q-factor from B&C
        self.avfpsq=float(DEFS["avfpsq"])
        self.aob=float(DEFS["aob"])  # alpha over beta  (Poisson solid)

        # high frequency magnitude correction
        self.HFEcorr=5.  # this is an old value from the legacy code. may need updating

        # processing tacers/energy
        self.cutoff=float(DEFS["cutoff"])  # factor by which to ignore data (values must be between mean/cutoff and mean*cutoff)
        self.snr=float(DEFS["snr"])  # minimum signal to noise ratio for processing
        self.pre_filt= (float(DEFS["pflp0"]), float(DEFS["pflp1"]), float(DEFS["pfhp1"]), float(DEFS["pfhp0"])) # low pass complete, low pass begin, high pass begin, high pass complete

        # Waveform source
        self.src = DEFS["src"]  # alternative 'IRIS'
    
class event:
    def __init__(self):
        #self.data="Existing"
        #self.data="New"  # 
        self.newData=True

        eloc = [38.8721,135.4642,367.70]
        etime= UTCDateTime(2021,9,29,8,37,6)   # North of Japan
        self.origin=[eloc,etime]
        self.ecount='00'
        self.iter='00'
        
        self.eventname=etime2name(etime,self.ecount)
        #additional event info
        phi,delta,lmbda=[70, 20, 90]
        self.focmech=[phi,delta,lmbda]

        #additional event info either determined or brought in
        self.Mi=0   # initial magnitude estimate regardless of type
        self.Mw=0
        self.Me=0
        self.ttime=0

def logmeanEnergy(E, Defaults=defaults(), **kwargs):
    """Calculate the log mean of energy E, excluding values below the cutoff"""
    cutoff = Defaults.cutoff
    
    logE = np.log10(E)
    last_vals = logE.iloc[-1]                # last value of each column
    
    # only 
    logEgt1 = logE.loc[:, last_vals > 1]
    last_valsgt1 = logEgt1.iloc[-1]

    mean_last_val = last_valsgt1.mean()         # mean of last values#logEclean = logE.loc[:, last_vals < mean_last_val+np.log10(cutoff) & last_vals > mean_last_val-np.log10(cutoff)]  # keep columns where last value < mean
    logEclean = logEgt1.loc[:, (last_valsgt1 < mean_last_val + np.log10(cutoff)) & (last_valsgt1 > mean_last_val - np.log10(cutoff))]
    print("keeping", logEclean.shape[1],"out of", logEgt1.shape[1],"traces due to cutoff=", cutoff)  # should be same as EBB.shape
    return 10 ** logEclean.mean(axis=1),logEclean.shape[1]  # mean of log values above cutoff

def bestWindow(E, windows,startTime=60, excludeLast=0, choice="MaxSlope", **kwargs):
    """
    Find the best window for growth or decay of the cumulative energy, E
    E = cumulative energy time series (pandas series)
    windows = array of window lengths to try (in seconds)
    startTime = time to start looking for windows (seconds)
    excludeLast = number of seconds at end of time series to exclude from consideration
    choice = "MaxSlope" or "MinMisfit" to choose best window using either the maximum slope or minimum misfit
    
    returns minMF_Window, minMF_t1, minMF_t2, minMF_slope, minMF_intercept, minMF, mresults
    """
    tlen=len(E)-excludeLast  # length of time series to use[seconds]. includes prePtime
    seq = np.arange(startTime, tlen + 1-excludeLast)
    # Fit growth of HF data
    mresults = []
    for window in windows: 
        results = []
        for t1 in seq[::1]:
            t2 = t1 + window    
            if t2 > tlen:  # keep it from running into a wall
                t2 = tlen
            wlen= t2 - t1  
            # Linear fit for EHFnormsum between t1 and t2
            x = np.arange(t1, t2)
            y = E.iloc[t1:t2]
            if len(x) >  2:  # Need at least 3 points for a fit
                coeffs = np.polyfit(x, y, 1)  # coeffs[0]=slope, coeffs[1]=intercept
                yfit = np.polyval(coeffs, x)
                misfit = np.sqrt(np.mean((y - yfit) ** 2))/(wlen-2)
            else:
                coeffs = [np.nan, np.nan]
                misfit = np.nan
            #print(f"t1={t1}, t2={t2}, slope={coeffs[0]:.4f}, intercept={coeffs[1]:.4f}, misfit={misfit:.4f}")
            results.append({
                "t1": t1,
                "t2": t2,
                "slope": coeffs[0],
                "intercept": coeffs[1],
                "misfit": misfit
            })
        if choice == "MinMisfit":
            idx = np.nanargmin([r["misfit"] for r in results])
        else:
        #elif choice == "MaxSlope":
            idx = np.nanargmax([r["slope"] for r in results])

        best_t1 = results[idx]["t1"]
        best_t2 = results[idx]["t2"]
        best_misfit = results[idx]["misfit"]
        best_slope = results[idx]["slope"]
        best_intercept = results[idx]["intercept"]
        #print(f"Peak slope: {peak_slope:.4f} at t1={peak_t1}, t2={peak_t2}, misfit={misfitPeak:.4e}, for window={window} seconds")
        mresults.append([window, best_t1, best_t2, best_slope, best_intercept, best_misfit])
    midx = np.nanargmin([m[5] for m in mresults])
    minMF_Window= mresults[midx][0]
    minMF_t1= mresults[midx][1]
    minMF_t2= mresults[midx][2]
    minMF_slope= mresults[midx][3]
    minMF_intercept= mresults[midx][4]
    minMF= mresults[midx][5]
    print(f"Best window: {minMF_Window} seconds, Peak slope: {minMF_slope:.2e} at t1={minMF_t1}, t2={minMF_t2}, misfit={minMF:.4e}")
    return minMF_Window, minMF_t1, minMF_t2, minMF_slope, minMF_intercept, minMF, mresults

def resultsWindow2Txo(r1,r2,prePtime):
    """
    r1 = results from the upslope fit (output of bestWindow)
    r2 = results from the downslope fit ("")
    prePtime=start time of waveforms relative to the theoretical P-wave (negative is before P)

    outputs Txo
    """
    min1Window, min1Peak_t1, min1Peak_t2, min1Peak_slope, min1Peak_intercept, min1MisfitPeak, m1results = r1
    min2Window, min2Peak_t1, min2Peak_t2, min2Peak_slope, min2Peak_intercept, min2MisfitPeak, m2results = r2
    return(min2Peak_intercept-min1Peak_intercept)/(min1Peak_slope-min2Peak_slope)+prePtime


def src2ergs(Defaults=defaults(), Event=event(), showPlots=False, **kwargs):
    """
    Run event processing from data retrieval to final results with plots saved in appropriate directories.  
    Heavy-lifting is performed by the classes (Defaults, and Event), so please be sure to locally define before running. 
    Those are where you define the basic data, processing, run and event parameters.

    Example:
    from rtergpy.run import defaults, event, etime2name, src2ergs
    from obspy import UTCDateTime
    
    Defaults=defaults()
    Event=event()

    Event.newData=False  # data may be already downloaded
    eloc = [-57.567,-25.032,47]  # event lat,lon,depth 
    etime= UTCDateTime(2021,8,12,18,32,52)  
    Event.eventname=etime2name(etime,ecount=Event.ecount)
    Event.origin=[eloc,etime]
    Event.focmech=[106, 26, 56] # phi,delta,lmbda

    Event.Mi  is initial magnitude and will be reset to Mw if it exists and is large

    Event.waveparams[1][1] #the amount of time to analyze following P-arrive will reset to larger values for big events
      # M8.0+  = 500 s
      # M8.5+ = 1000 s

    src2ergs(Defaults=Defaults,Event=Event)
    """
    from rtergpy.waveforms import getwaves,ErgsFromWaves,loadwaves,gmeanCut,tacer,tacerstats
    from rtergpy.waveforms import trstat2pd,e2Me,eventdir, load_seismoGNSS_waves
    from rtergpy.plotting import tacerplot,Edistplot,Ehistogram,Eazplot, stationEmapPygmt, fbandlabels, Efluxplots,ETxoplot
    from scipy.interpolate import interp1d
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import os
    
    eventname=Event.eventname
    
    # reset some parameters based on EQ size
    # set initial magnitude to Mw if Mw is big and Mi is not
    if (Event.Mi < 1) and (Event.Mw > 1) :
        Event.Mi = Event.Mw
    
    # extend duration of data to get and run on for the largest earthquakes
    if Event.Mi > 8.5 : 
        Defaults.waveparams[1][1]=1000  # very bit event so go longer
    elif Event.Mi > 8 : 
        Defaults.waveparams[1][1]=500

    if Event.focmech is False :
        Event.focmech = [0,30,90] #using generic mechanism (thrust)
        print("WARNING: No mechanism found, setting to a default")
    
    if getattr(Event, "is_seismogeodetic", False):
        print("Using seismo-geodetic waveforms provided externally (no download).")
        st, df = load_seismoGNSS_waves(Defaults=Defaults, Event=Event, event_mseed=Event.mseed, df_path=Event.dfpath)
    
        if st is None or df is None:
            raise ValueError("Seismogeodetic mode requires st and df to be passed to src2ergs.")
    else:
        if Event.newData:
            print("Getting waveforms")
            st,df=getwaves(Defaults=Defaults,Event=Event)
        else:
            print("Loading locally stored waveforms")
            try:
                st,df=loadwaves(Defaults=Defaults,Event=Event) 
            except:
                print("Couldn't load data for "+eventname+". Attempting to download:")
                st,df=getwaves(Defaults=Defaults,Event=Event)
        if len(st) == 0:
            raise Exception("ERROR: No waveforms retreived.") 
    #runiter= iterate()
        # runiter.count gives current count
        # runiter.str() gives as a 2digit string leading zero
        # runiter.step() increases iteration by 1

# remove any additional locations at same site (or site repeats!)
    prePtime=Defaults.waveparams[1][0]
    STATCHAN0=''
    for tr in st:
        STATCHAN=str(tr.stats.network)+str(tr.stats.station)
        trWindow=UTCDateTime(tr.stats.endtime)-UTCDateTime(tr.stats.starttime)
        nDataExpected=trWindow*tr.stats.sampling_rate+1

        if STATCHAN == STATCHAN0:  # multiple sensors at same location
            st.remove(tr)
        elif trWindow <  0-prePtime:  # too short a window (prePtime is negative)
            st.remove(tr)
        elif tr.stats.npts < nDataExpected:  # less data than expected
            st.remove(tr)
        STATCHAN0=STATCHAN

    # work in event directory 
    edirit,origwd=eventdir(Defaults=Defaults,Event=Event,create=False,cd=True)
    df.eventdir = edirit
    # create an array of station data 
    trdf=pd.DataFrame()
    trdf = pd.concat([trstat2pd(tr) for tr in st], ignore_index=True)
    #for tr in st:
    #    trdf=trdf.append(trstat2pd(tr),ignore_index=True)

    print("Calculating Energy growth with time ")
    # Calculate energies over all waves
    EBB,EHF,Emd=ErgsFromWaves(st,Defaults,Event)   
    print("Length EBB and EHF", len(EBB),len(EHF))
    # get tacer values and fist derivatives 
    print("Calculating TACER Values")
    #kern=10  # 1/2 width of gauss @ 1sig
    kern=Defaults.smoothkern
    EBBSmooth=EBB.rolling(kern,win_type='gaussian',center=True, closed='both').mean(std=kern/2)
    dEBBdt=EBB.diff()
    dEBBdtSmooth=EBBSmooth.diff()

    EHFSmooth=EHF.rolling(kern,win_type='gaussian',center=True, closed='both').mean(std=kern/2)
    dEHFdt=EHF.diff()
    dEHFdtSmooth=EHFSmooth.diff()

    tacerBB=tacer(dEBBdtSmooth,prePtime=prePtime)
    tacerHF=tacer(dEHFdtSmooth,prePtime=prePtime)
    ttimes,meds = tacerstats(tacerHF) 
    ttimeHF,ttimeHF25,ttimeHF75=ttimes
    print("Median Tacer time = %.1f -/+ %.1f/%.1f s (25/75th percentile)" %(ttimeHF,ttimeHF25,ttimeHF75))
    # === Minimal station alignment (no architectural changes) ===
    # Source of truth = columns present in EBB (i.e., stations that survived energy)
    station_cols = list(EBB.columns)

    # Reindex EHF to EBB columns (usually already matches; no-op for BB)
    EHF = EHF.reindex(columns=station_cols)

    # Reindex metadata tables to the same stations and order
    if isinstance(Emd, pd.DataFrame) and 'netstatchan' in Emd.columns:
        Emd = (
            Emd.drop_duplicates(subset=['netstatchan'], keep='first')
            .set_index('netstatchan')
            .reindex(station_cols)
            .reset_index()
        )

    if isinstance(trdf, pd.DataFrame) and 'netstatchan' in trdf.columns:
        trdf = (
            trdf.drop_duplicates(subset=['netstatchan'], keep='first')
                .set_index('netstatchan')
                .reindex(station_cols)
                .reset_index()
        )

    # If use per-station meds indexed by the original column order, trim to the same stations
    #    (assumes meds rows correspond to EBB/EHF columns in order)
    if isinstance(meds, pd.DataFrame) and meds.shape[0] != len(station_cols):
        if 'netstatchan' in meds.columns:
            mask = meds['netstatchan'].isin(station_cols)
            meds = meds.loc[mask].reset_index(drop=True)
        else:
            meds = meds.iloc[:len(station_cols)].reset_index(drop=True)

    intval=int(ttimeHF-prePtime)  # median tacer time from start of waveforms
    ebbmedtac=np.array(EBB.iloc[intval])  # list of EBB values at *median* tacer
    ehfmedtac=np.array(EHF.iloc[intval])  # list of EHF values at *median* tacer
    
    # corrected values for focmech
    ebbcorrmedtac=np.array(list(ebbmedtac*np.array(Emd.est2corr)))
    ehfcorrmedtac=np.array(list(ehfmedtac*np.array(Emd.est2corr)))
    
    # Energy values at per station tacer
    ebbpertac=[]
    ehfpertac=[]
    for i in range(0,len(meds)):
        ebbpertac = np.append(ebbpertac, EBB.iloc[meds['time at max'][i]-prePtime][i])
        ehfpertac = np.append(ehfpertac, EHF.iloc[meds['time at max'][i]-prePtime][i])

    ebbcorrpertac=np.array(list(ebbpertac*np.array(Emd.est2corr)))
    ehfcorrpertac=np.array(list(ehfpertac*np.array(Emd.est2corr)))


    # cutoff=15 # 15x +/-  # moved into defaults
    cutoff=Defaults.cutoff
    labelbb,labelhf=fbandlabels([Emd.fbandBB[0],Emd.fbandHF[0]])

    print("From Median Tacer: --------------------------")
    ebbmedtacmean,keepbb=gmeanCut(ebbmedtac,cutoff=cutoff)
    ehfmedtacmean,keephf=gmeanCut(ehfmedtac,cutoff=cutoff)
    ebbmedtacmeanerr10=np.std(np.log10(keepbb))
    ehfmedtacmeanerr10=np.std(np.log10(keephf))
    print("  Mean BB Energy (Estimated)= %.2e [Me %.2f]" %(ebbmedtacmean,e2Me(ebbmedtacmean)))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelbb,len(keepbb), ebbmedtacmean, ebbmedtacmeanerr10))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelhf,len(keephf), ehfmedtacmean, ehfmedtacmeanerr10))
    ebbcorrmedtacmean,keepbb=gmeanCut(ebbcorrmedtac,cutoff=cutoff)
    ehfcorrmedtacmean,keephf=gmeanCut(ehfcorrmedtac,cutoff=cutoff)
    ebbcorrmedtacmeanerr10=np.std(np.log10(keepbb))
    ehfcorrmedtacmeanerr10=np.std(np.log10(keephf))
    print("  Mean BB Energy (FM corrected) = %.2e [Me %.2f]" %(ebbcorrmedtacmean,e2Me(ebbcorrmedtacmean)))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelbb,len(keepbb), ebbcorrmedtacmean, ebbcorrmedtacmeanerr10))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelhf,len(keephf), ehfcorrmedtacmean, ehfcorrmedtacmeanerr10))

    print("From Per-Station Tacer: ---------------------")
    ebbpertacmean,keepbb=gmeanCut(ebbpertac,cutoff=cutoff)
    ehfpertacmean,keephf=gmeanCut(ehfpertac,cutoff=cutoff)
    ebbpertacmeanerr10=np.std(np.log10(keepbb))
    ehfpertacmeanerr10=np.std(np.log10(keephf))
    print("  Mean BB Energy (Estimated)= %.2e [Me %.2f]" %(ebbpertacmean,e2Me(ebbpertacmean)))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelbb,len(keepbb), ebbpertacmean, ebbpertacmeanerr10))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelhf,len(keephf), ehfpertacmean, ehfpertacmeanerr10))
    ebbcorrpertacmean,keepbb=gmeanCut(ebbcorrpertac,cutoff=cutoff)
    ehfcorrpertacmean,keephf=gmeanCut(ehfcorrpertac,cutoff=cutoff)
    ebbcorrpertacmeanerr10=np.std(np.log10(keepbb))
    ehfcorrpertacmeanerr10=np.std(np.log10(keephf))
    print("  Mean BB Energy (FM corrected) = %.2e [Me %.2f]" %(ebbcorrpertacmean,e2Me(ebbcorrpertacmean)))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelbb,len(keepbb), ebbcorrpertacmean, ebbcorrpertacmeanerr10))
    print("    %s: %d traces, %.2e +- 10^%.2f [J]" %(labelhf,len(keephf), ehfcorrpertacmean, ehfcorrpertacmeanerr10))
    #print("Mean BB Energy (FM corrected) = %.2e [Me %.2f]" %(emeanbbcorr,e2Me(emeanbbcorr)))

    # Saving information  ################################
    # create dataframe with Event based results
    results=pd.DataFrame({
        "eventname":eventname,"iteration":Event.iter,   # name/run
        "etime":Event.origin[1],"elat":Event.origin[0][0],"elon":Event.origin[0][1],"edepth":Event.origin[0][2], "focmech":[Event.focmech],   # Event
        "network":Defaults.network, "chan":Defaults.chan, "stationrange":[Defaults.stationrange], "nstats":len(trdf), #stations 
        "fbands":[Defaults.waveparams[0]], "pwindow":[Defaults.waveparams[1]],  # wave params
        "eventdir":df.eventdir, "modtime":UTCDateTime(),      # where and when processed
        "cutoff":Defaults.cutoff, "ttimes":[ttimes],  # Results (and below)
        "ebbmedtacmean":ebbmedtacmean, "STD10(medtac)":ebbmedtacmeanerr10, "Me(medtac)":e2Me(ebbmedtacmean), "ehfmedtacmean":ehfmedtacmean, "STD10(hfmedtac)":ehfmedtacmeanerr10,
        "ebbcorrmedtacmean":ebbcorrmedtacmean, "STD10(corr)":ebbcorrmedtacmeanerr10, "Me(corr)":e2Me(ebbcorrmedtacmean), "ehfcorrmedtacmean":ehfcorrmedtacmean, "STD10(hfcorr)":ehfcorrmedtacmeanerr10,
        "ebbpertacmean":ebbpertacmean, "STD10(per)":ebbpertacmeanerr10, "Me(per)":e2Me(ebbpertacmean), "ehfpertacmean":ehfpertacmean, "STD10(hfper)":ehfpertacmeanerr10,
        "ebbcorrpertacmean":ebbcorrpertacmean, "STD10(percorr)":ebbcorrpertacmeanerr10, "Me(corrper)":e2Me(ebbcorrpertacmean), "ehfcorrpertacmean":ehfcorrpertacmean, "STD10(hfpercorr)":ehfcorrpertacmeanerr10
        }, dtype=object)
    
    # time-series energy results
    Etimeseries=pd.concat([EBB,EHF,EBBSmooth,EHFSmooth,dEBBdtSmooth,dEHFdtSmooth,tacerBB,tacerHF],
        keys=["EBB","EHF","EBBSmooth","EHFSmooth","dEBBdtSmooth","dEHFdtSmooth","tacerBB","tacerHF"])
        # individual key can be extracted using (e.g. Energy.loc["EHF"])

    # per-station information
    ETace=pd.DataFrame({'tacer':meds['time at max'],
        'ebbmedtac':ebbmedtac, 'ehfmedtac':ehfmedtac,
        'ebbcorrmedtac':ebbcorrmedtac, 'ehfcorrmedtac':ehfcorrmedtac,
        'ebbpertac':ebbpertac, 'ehfpertac':ehfpertac,
        'ebbcorrpertac':ebbcorrpertac, 'ehfcorrpertac':ehfcorrpertac
        })
    ETace=ETace.reset_index(drop=True)
    StationTacer=pd.concat([trdf,Emd[["estFgP2","FgP2","est2corr"]],ETace],axis=1)

    Event.Me=e2Me(ebbpertacmean)
    Event.ttime=ttimeHF

    #.  Create Legecy Txo, Ebb and Ehf from Txo
    try: 
        HFEcorr=Defaults.HFEcorr
        EBBTxo = 0
        EHFTxo = 0
        Txo = 0

        # should create results that are the same/similar to before.  We would average only results withing 15x the original geometric mean.
        print("Calculating log mean of EBB and EHF energies for cumulative growth curves")
        EBBlmean,nBBlmean=logmeanEnergy(EBB)
        EHFlmean,nHFlmean=logmeanEnergy(EHF)
        Event.nBBlmean=nBBlmean
        Event.nHFlmean=nHFlmean

        if (nBBlmean > 0) and (nHFlmean > 0) : 
            # Extract best window fit info
            fullwindow=int(Defaults.waveparams[1][1]+Defaults.waveparams[1][0])
            halfwindow=int(fullwindow/2)
            quarterwindow=int(fullwindow/4)
            
            windowUp= np.arange(3, quarterwindow, 1) # use at least 3 second and up to 1/4 of the total time series for fitting
            upResults=bestWindow(EHFlmean, windows=windowUp, startTime=0-prePtime, excludeLast=halfwindow, choice="MaxSlope",  minwindow=3) # find the best fit and window for the up-slope (controlled in part by prePtime and window choices)
            startTimeDown=upResults[2]+upResults[0] # end of upslope fit + window used 
            windowDown= np.arange(halfwindow,fullwindow,5) # tries to use as much data as possible after the upslope
            downResults=bestWindow(EHFlmean, windows=windowDown, startTime=startTimeDown, excludeLast=60, choice="MinMisfit") # find the best fit and window for the down-slope (controlled in part by prePtime and window choices)

            # Extract best window fit info
            Txo=resultsWindow2Txo(upResults,downResults,prePtime)
            print(f"Txo = {Txo:.2f} seconds (crossover time)")

            # Energy values at Txo
            x=np.arange(0, len(EHFlmean))
            EBBlmean_func = interp1d(x + prePtime, EBBlmean, bounds_error=False, fill_value="extrapolate")
            EBBTxo = EBBlmean_func(Txo)
            print(f"EBBlmean at Txo ({Txo:.2f} s): {EBBTxo:.2e} J, (MeBB {e2Me(EBBTxo):.2f})")

            EHFlmean_func = interp1d(x + prePtime, EHFlmean, bounds_error=False, fill_value="extrapolate")
            EHFTxo = EHFlmean_func(Txo)
            print(f"EHFlmean at Txo ({Txo:.2f} s): {EHFTxo:.2e} J, (MeHF {e2Me(EHFTxo,eCorrection=HFEcorr):.2f})")
        else:
            print("Notice: Insufficient data kept for logmean cumulative energy growth curves. nBBlmean = %i, nHFlmean = %" % ( nBBlmean,nHFlmean))
    except: 
        print("Notice: Could not create Legacy Txo, EBBTxo, or EHFTxo")

    Event.EBBTxo=EBBTxo
    Event.EHFTxo=EHFTxo
    Event.Txo=Txo

    # Create plots  
    print("Making figures\n")

    if not os.path.exists('figs'):   # create and go into pkls dir
        os.mkdir('figs')
    
    os.chdir('figs')
    if not showPlots:
        mpl.use('Agg')  # needed to plot without using the X-session
    # individual plot runs    
    droppcts=[0.5,.25,0.1]
    try: 
        eloc=Event.origin[0]
        st.plot(type='section', dist_degree=True, ev_coord=[eloc[0], eloc[1]], orientation='horizontal', scale=3, outfile=f"Moveout_{Event.eventname}.png");
    except:
            print("ERROR: Moveout Plot for "+eventname+":",e)
    try:
        droptimes=Efluxplots(dEHFdtSmooth, trdf, Event=Event, Defaults=Defaults, pcts=droppcts, show=showPlots)    
    except:
        print("ERROR: Eflux plots for "+eventname+":",e)
    try:
        tacerplot(tacerHF,trdf,ttimes,meds,eventname,show=showPlots)
    except:
        print("ERROR: Tacer plot for  "+eventname+":",e)
    try:    
        Edistplot(EBB,EHF,Emd,trdf,eventname,ttimeHF, prePtime=prePtime,show=showPlots,cutoff=cutoff)
    except:
        print("ERROR: Edistance plot for  "+eventname+":",e)
    try:
        Eazplot(EBB,EHF,Emd,trdf,eventname,ttimeHF, prePtime=prePtime,show=showPlots,cutoff=cutoff)
    except:
        print("ERROR: E-AZ plot for  "+eventname+":",e)
    try:
        Ehistogram(EBB,EHF,Emd,eventname,ttimeHF, prePtime=prePtime,show=showPlots,cutoff=cutoff)
    except:
        print("ERROR: E histogram for  "+eventname+":",e)
    try:
        if Txo > 0 :
            ETxoplot(EBBlmean,EHFlmean,upResults,downResults,Event=Event,Defaults=Defaults,show=False)
        else:
            print("Notice: Skipping ETxo plotting, Txo = 0")
    except:
        print("ERROR: E Txo plot for  "+eventname+":",e)

    mpl.pyplot.close('all')  # they don't close themselves

    try:
        stationEmapPygmt(EBB,trdf,ttimeHF,Event=Event,Defaults=Defaults,prefix="BB_",show=showPlots)
        stationEmapPygmt(EHF,trdf,ttimeHF,Event=Event,Defaults=Defaults,prefix="HF_",show=showPlots)
    #    stationEmapPygmt(EBB,Event.origin[0],trdf,eventname,ttimeHF, prePtime=prePtime,cutoff=15,itername=Event.iter,show=showPlots)
    except IOError as err:
        print("ERROR: Map Plot for "+eventname+":",e)
    os.chdir(edirit)

    # Saving information  ################################
    # create dataframe with Event based results
    results=pd.DataFrame({
        "eventname":eventname,"iteration":Event.iter,   # name/run
        "etime":Event.origin[1],"elat":Event.origin[0][0],"elon":Event.origin[0][1],"edepth":Event.origin[0][2], "focmech":[Event.focmech],   # Event
        "network":Defaults.network, "chan":Defaults.chan, "stationrange":[Defaults.stationrange], "nstats":len(trdf), #stations 
        "fbands":[Defaults.waveparams[0]], "pwindow":[Defaults.waveparams[1]],  # wave params
        "eventdir":df.eventdir, "modtime":UTCDateTime(),      # where and when processed
        "cutoff":Defaults.cutoff, "ttimes":[ttimes],  # Results (and below)
        "ebbmedtacmean":ebbmedtacmean, "STD10(medtac)":ebbmedtacmeanerr10, "Me(medtac)":e2Me(ebbmedtacmean), 
        "ehfmedtacmean":ehfmedtacmean, "STD10(hfmedtac)":ehfmedtacmeanerr10,
        "ebbcorrmedtacmean":ebbcorrmedtacmean, "STD10(corr)":ebbcorrmedtacmeanerr10, "Me(corr)":e2Me(ebbcorrmedtacmean), 
        "ehfcorrmedtacmean":ehfcorrmedtacmean, "STD10(hfcorr)":ehfcorrmedtacmeanerr10,
        "ebbpertacmean":ebbpertacmean, "STD10(per)":ebbpertacmeanerr10, "Me(per)":Event.Me, 
        "ehfpertacmean":ehfpertacmean, "STD10(hfper)":ehfpertacmeanerr10,
        "ebbcorrpertacmean":ebbcorrpertacmean, "STD10(percorr)":ebbcorrpertacmeanerr10, "Me(corrper)":e2Me(ebbcorrpertacmean), 
        "ehfcorrpertacmean":ehfcorrpertacmean, "STD10(hfpercorr)":ehfcorrpertacmeanerr10,
        "Ebb(Txo)": Event.EBBTxo, "Me(Txo)":e2Me(Event.EBBTxo), "EHF(Txo)":Event.EHFTxo, "Mehf(Txo)":e2Me(Event.EHFTxo), "Txo":Event.Txo,
        "Droptimes": [droptimes]
        }, dtype=object)

    # time-series energy results for each station
    Etimeseries=pd.concat([EBB,EHF,EBBSmooth,EHFSmooth,dEBBdtSmooth,dEHFdtSmooth,tacerBB,tacerHF],
        keys=["EBB","EHF","EBBSmooth","EHFSmooth","dEBBdtSmooth","dEHFdtSmooth","tacerBB","tacerHF"])
        # individual key can be extracted using (e.g. Energy.loc["EHF"])

    # Cumulative Energy Time Series using the logrithmic mean 
    Ecumtimeseries=pd.concat([EBBlmean,EHFlmean], keys=["EBBlmean","EHFlmean"], axis=1)

    # per-station information
    ETace=pd.DataFrame({'tacer':meds['time at max'],
        'ebbmedtac':ebbmedtac, 'ehfmedtac':ehfmedtac,
        'ebbcorrmedtac':ebbcorrmedtac, 'ehfcorrmedtac':ehfcorrmedtac,
        'ebbpertac':ebbpertac, 'ehfpertac':ehfpertac,
        'ebbcorrpertac':ebbcorrpertac, 'ehfcorrpertac':ehfcorrpertac
        })
    ETace=ETace.reset_index(drop=True)
    StationTacer=pd.concat([trdf,Emd[["estFgP2","FgP2","est2corr"]],ETace],axis=1)

    # save results to files
    try:
        print("writing results\n")
        # csv  
        if not os.path.exists('csvs'):   # create and go into pkls dir
            os.mkdir('csvs')

        os.chdir('csvs')  
        results.to_csv("Results_"+eventname+".csv")
        Etimeseries.to_csv("EStationTimeSeries_"+eventname+".csv")
        StationTacer.to_csv("ETacer_"+eventname+".csv")
        Ecumtimeseries.to_csv("ECumulativeTimeSeries_"+eventname+".csv") 
        os.chdir(edirit)
        # pkls 
        if not os.path.exists('pkls'):   # create and go into pkls dir
            os.mkdir('pkls')
        os.chdir('pkls')
        results.to_pickle("Results_"+eventname+".pkl")
        Etimeseries.to_pickle("EStationTimeSeries_"+eventname+".pkl")
        StationTacer.to_pickle("ETacer_"+eventname+".pkl")
        Ecumtimeseries.to_csv("ECumulativeTimeSeries_"+eventname+".csv") 
        os.chdir(edirit)
    except:
        print("ERROR: writing results for"+eventname)
    os.chdir(origwd)  # go back to old directory

def mergeResults(Defaults=defaults(), iteration='00', **kwargs):
    """
    Reads all processed event information and returns a master dataframe of summary result information
    """ 
    import glob
    import pandas as pd
    
    
    files= glob.glob(Defaults.edirbase +'/[12]???/[12]*/'+iteration+'/pkls/Results*.pkl')
    prior='' 
    for file in files:
        if not prior: # create first time
            df=pd.read_pickle(file)
            prior=1  # no longer first
        else:
            dflocal=pd.read_pickle(file)
            df=df.append(dflocal,ignore_index=True)
    # replace list ttimes wih columsn of individual values
    dfttimes=pd.DataFrame(df['ttimes'].to_list(), columns = ['tacer', 't25', 't75'])
    df.drop('ttimes', axis=1, inplace=True)
    df.insert(16, 'tacer', dfttimes['tacer'],True)
    df.insert(17, 't25', dfttimes['t25'],True)
    df.insert(18, 't75', dfttimes['t75'],True)
    df.sort_values(by=['eventname'],inplace=True,ignore_index=True)  # results should be time sorted now
    return df
