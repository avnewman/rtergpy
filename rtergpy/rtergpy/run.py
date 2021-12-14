#from rtergpy.waveforms import etime2name
#from rtergpy.waveforms import tacerstats,trstat2pd,e2Me,process_waves,eventdir,iterate
#from rtergpy.plotting import tacerplot,Edistplot,Ehistogram,Eazplot,stationEmapBasemap, stationEmapPygmt
from obspy import UTCDateTime
from tqdm import tqdm
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
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
        
        # frequency range
        f1min=float(DEFS["f1min"]) ; f1max=float(DEFS["f1max"]) # BB
        f2min=float(DEFS["f2min"]) ; f2max=float(DEFS["f2max"]) # HF 
        fbands=[[f1min,f1max],[f2min,f2max]]
        tstep=int(DEFS["tstep"]) # increment
        prePtime=int(DEFS["prePtime"])
        postPtime=int(DEFS["postPtime"])
        pwindow=[prePtime,postPtime]
        self.waveparams=[fbands,pwindow,tstep]
        self.resample=int(DEFS["resample"])  # samples per second
        self.smoothkern=int(DEFS["smoothkern"])  # kernel for gaussian smoothing of data for duration estimates (1/2 on each side)

        # earth params
        self.rearth=float(DEFS["rearth"]) # in meters
        self.qbc=float(DEFS["qbc"]) #q-factor from B&C
        self.avfpsq=float(DEFS["avfpsq"])
        self.aob=float(DEFS["aob"])  # alpha over beta  (Poisson solid)

        # processing tacers/energy
        self.cutoff=float(DEFS["cutoff"])  # factor by which to ignore data (values must be between mean/cutoff and mean*cutoff)

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
        self.Mw=0
        self.Me=0
        self.ttime=0


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

    src2ergs(Defaults=Defaults,Event=Event)
    """
    from rtergpy.waveforms import getwaves,ErgsFromWaves,loadwaves,gmeanCut,tacer,tacerstats
    from rtergpy.waveforms import trstat2pd,e2Me,eventdir
    from rtergpy.plotting import tacerplot,Edistplot,Ehistogram,Eazplot, stationEmapPygmt, fbandlabels, Efluxplots
    import numpy as np
    import pandas as pd
    import matplotlib as mpl
    import os
    
    eventname=Event.eventname

    if Event.newData:
        print("Getting waveforms")
        #st,df=getwaves(eloc,etime,pwindow=Def.waveparams[1],rads=Def.stationrange)
        st,df=getwaves(Defaults=Defaults,Event=Event)
    else:
        print("Loading locally stored waveforms")
        try:
            st,df=loadwaves(Defaults=Defaults,Event=Event) 
        except:
            print("Couldn't load data for "+eventname+". Quitting.")
            exit(1)
    if len(st) == 0:
        raise Exception("ERROR: No waveforms retreived.") 
    #runiter= iterate()
        # runiter.count gives current count
        # runiter.str() gives as a 2digit string leading zero
        # runiter.step() increases iteration by 1

# remove any additional locations at same site (or site repeats!)
    STATCHAN0=''
    for tr in st:
        STATCHAN=str(tr.stats.network)+str(tr.stats.station)
        if STATCHAN == STATCHAN0:
            st.remove(tr)
        STATCHAN0=STATCHAN

    # work in event directory 
    edirit,origwd=eventdir(Defaults=Defaults,Event=Event,create=False,cd=True)

    # create an array of station data 
    trdf=pd.DataFrame()
    for tr in st:
        trdf=trdf.append(trstat2pd(tr),ignore_index=True)

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

    prePtime=Defaults.waveparams[1][0]
    tacerBB=tacer(dEBBdtSmooth,prePtime=prePtime)
    tacerHF=tacer(dEHFdtSmooth,prePtime=prePtime)
    ttimes,meds = tacerstats(tacerHF) 
    ttimeHF,ttimeHF25,ttimeHF75=ttimes
    print("Median Tacer time = %.1f -/+ %.1f/%.1f s (25/75th percentile)" %(ttimeHF,ttimeHF25,ttimeHF75))
    
    intval=int(ttimeHF-prePtime)  # median tacer time from start of waveforms
    ebbmedtac=np.array(EBB.iloc[intval])  # list of EBB values at *median* tacer
    ehfmedtac=np.array(EHF.iloc[intval])  # list of EBB values at *median* tacer
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
    labelbb,labelhf=fbandlabels(Emd)

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

    # Create plots  
    try:
        print("Making figures\n")
        if not os.path.exists('figs'):   # create and go into pkls dir
            os.mkdir('figs')
        os.chdir('figs')
        if not showPlots:
            mpl.use('Agg')  # needed to plot without using the X-session
        # individual plot runs    
        droppcts=[0.5,.25,0.1]
        droptimes=Efluxplots(dEHFdtSmooth, trdf, Event=Event, Defaults=Defaults, pcts=droppcts, show=showPlots)    
        tacerplot(tacerHF,trdf,ttimes,meds,eventname,show=showPlots)
        Edistplot(EBB,EHF,Emd,trdf,eventname,ttimeHF, prePtime=prePtime,show=showPlots,cutoff=cutoff)
        Eazplot(EBB,EHF,Emd,trdf,eventname,ttimeHF, prePtime=prePtime,show=showPlots,cutoff=cutoff)
        Ehistogram(EBB,EHF,Emd,eventname,ttimeHF, prePtime=prePtime,show=showPlots,cutoff=cutoff)
        stationEmapPygmt(EBB,Event.origin[0],trdf,eventname,ttimeHF, prePtime=prePtime,cutoff=15,itername=Event.iter,show=showPlots)
        os.chdir('..')
        mpl.pyplot.close('all')  # they don't close themselves
    except:
        print("ERROR: plotting results for "+eventname)

    results["Droptimes"]=[droptimes]

    # save results to files
    try:
        print("writing results\n")
        # csv    
        results.to_csv("Results_"+eventname+".csv")
        Etimeseries.to_csv("Etimeseries_"+eventname+".csv")
        StationTacer.to_csv("ETacer_"+eventname+".csv") 
        # pkls 
        if not os.path.exists('pkls'):   # create and go into pkls dir
            os.mkdir('pkls')
        os.chdir('pkls')
        results.to_pickle("Results_"+eventname+".pkl")
        Etimeseries.to_pickle("Etimeseries_"+eventname+".pkl")
        StationTacer.to_pickle("Etacer_"+eventname+".pkl")
        os.chdir('..')
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
