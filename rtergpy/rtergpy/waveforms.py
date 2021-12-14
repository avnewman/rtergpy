from obspy.core.stream import Stream
from obspy.clients.neic import Client as nClient
from obspy.clients.fdsn import Client as fdsnClient
from obspy import UTCDateTime
from obspy.geodetics.base import locations2degrees as l2d
from obspy.geodetics.base import degrees2kilometers as d2km
from obspy.geodetics.base import kilometers2degrees as km2d
from obspy.geodetics.base import gps2dist_azimuth as ll2az
from obspy.taup import TauPyModel
from numpy import sin,cos,arcsin,sqrt,abs,pi,log10,exp
from scipy.fftpack import fft,ifft
from scipy.stats import gmean 
from tqdm import tqdm
from compress_pickle import dump as cpkldump # reading/writing compressed pickles
#from compress_pickle import load as cpklload # reading/writing compressed pickles
import os
import numpy as np
import pandas as pd
model = TauPyModel(model="iasp91")
from rtergpy.run import defaults, event
Defaults=defaults()
Event=event()
#edirbase="/Users/anewman/Documents/Projects/rterg/processing/Examples/events"
#network = "CU,GT,IC,II,IU" # GSN network
#stationrange=[25,80]  # distance in degrees
#chan="BHZ"
#prePtime=-60; postPtime=300
#pwindow=[prePtime,postPtime]
#resample=10  # samples per second

def process_waves(st,taper=0.05, freqmin=0.01, freqmax=2, **kwargs):
    """
    performs default waveform processing, on a copy, that is returned.
    """
    # waveform processing
    stp=st.copy()  # create backup
    # process data
    stp.detrend(type='polynomial', order=5) # pre-instrument removal
    stp.taper(taper)
    stp.remove_response(output="VEL")
    stp.detrend(type='polynomial', order=5) # post-instrument removal
    stp.taper(taper)

    if freqmin < 0.00001:  # esentially 0
        stp.filter("lowpass", freq=freqmax)
    else:
        stp.filter("bandpass", freqmin=freqmin, freqmax=freqmax)
    return stp
    
def theorPinfo(eloc,etime,sloc):
    """
    Algorithm to use taup event, and station location info to calculate theoretical P-time
    and other ray info.  Produces
    Ptime in UTC
    Ptakeoff angle in degrees
    Prayparameter in sec/rads
    Pincidence angle in degrees
    Distance in degrees
    A.V. Newman Mon Jul 26 15:26:35 EDT 2021
    """
    elat,elon,edep=eloc
    slat,slon,sheight=sloc
    arrivals=model.get_travel_times_geo(
        source_depth_in_km=edep, source_latitude_in_deg=elat, source_longitude_in_deg=elon,
        receiver_latitude_in_deg=slat, receiver_longitude_in_deg=slon,
        phase_list="P")

    ptime=etime + arrivals[0].time
    toa=arrivals[0].takeoff_angle
    rayp=arrivals[0].ray_param
    inc=arrivals[0].incident_angle
    distdeg=arrivals[0].distance
    return  ptime,toa,rayp,inc,distdeg


def get_respinv(network,eloc,etime,rads,chan):
    """
    Algorithm to use Obspy's metadata to pull response and other metadata.  Returns an apporpiate
    inventory class.
    A.V. Newman Mon Jul 26 15:26:35 EDT 2021
    """
    fclient = fdsnClient()

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

def downloadwaves(inventory,eloc,etime,pwindow=Defaults.waveparams[1],src=Defaults.src, **kwargs):
    """
    Module to use Obspy's metadata from get_respinv (fsdn) to pull data from
    the NEIC server around the P-wave theoretical time.  
    Returns an streams with added station location, event distance metadata.

    A.V. Newman Mon Jul 26 15:26:35 EDT 2021
    """
    
    if src == 'NEIC':
        print("Getting waves from NEIC")
        nclient=nClient()
    #elif src == 'ISC':
    #    print("Getting waves from ISC")
    #    fclient=fdsnClient("ISC")
    else:
        print("Getting waves from IRIS")
        fclient=fdsnClient("IRIS")

    st=Stream()

    # run on all channels in inventory
    for chan in tqdm(inventory.get_contents().get("channels")):
        slat,slon,sz,sldep=inventory.get_coordinates(chan).values()
        sheight=sz-sldep
        sloc=slat,slon,sheight
        Ptime,Ptoa,Prayp,Pinc,distdeg=theorPinfo(eloc,etime,sloc)
        distmeters,az,baz=ll2az(eloc[0],eloc[1],slat,slon)
        neti,stati,loci,chani=chan.split(".")
        stlocal='' # start with empty field in case first wave fails
        try:
            if src == 'NEIC':
                stlocal=nclient.get_waveforms(neti,stati,loci,chani,Ptime+pwindow[0],Ptime+pwindow[1])
            else:
                stlocal=fclient.get_waveforms(neti,stati,loci,chani,Ptime+pwindow[0],Ptime+pwindow[1], minimumlength=120, longestonly=True)
            #print("%s.%s.%s.%s downloaded. Continuing.." %(neti,stati,loci,chani))
        except:
            print("%s.%s.%s.%s failed to download. Continuing.." %(neti,stati,loci,chani))
        if stlocal:
            # add station coordinates
            try:
                stlocal[0].stats.coordinates= {'latitude': slat, 'longitude': slon}
                stlocal[0].stats.distance=distmeters;  # distance should be reported in meters!
                stlocal[0].stats.phasePtime=Ptime;  # UTC time of P-arrival
                stlocal[0].stats.ptoa=Ptoa;  # take-off
                stlocal[0].stats.prayp=Prayp;  # ray parameter
                stlocal[0].stats.pinc=Pinc;  # ray parameter
                stlocal[0].stats.distdeg=distdeg;  # UTC time of P-arrival
                stlocal[0].stats.az=az;  # azimuth
                stlocal[0].stats.baz=baz;  # back-azimuth
                st+=stlocal[0]
            except:
                print ("Channel ", chan, " not added...missing metadata")
        stlocal='' # clear field 
        #else:
            #print("- skip - :", chan)
        
    st.attach_response(inventory) # include response information
    # finally remove any tr that doesn't include response info
    for tr in st:
        try: 
            tr.stats.response
        except:
            print("No Response info for below. removing.\n", tr)
            st.remove(tr)
    return st 

def tstar(f):
  """
  Tstar as a function of frequency for teleseismic shallow EQs
  Digitized from Choy and Boatwright 1995 for shallow events at teleseismic distances
  (which they derived from Choy and Cormier (1986) -which shows only a 600km depth event)
  """
  if (f < 0.1) :
    tst=0.9-0.1*log10(f)
  elif (f < 1) :
    tst=0.5-0.5*log10(f)
  else :
    tst=0.5-.10*log10(f)  # TODO check values
  return tst

def gttP(depkm,distdeg):
  """
  returns travel time in seconds given event depth (km) and distance (degrees)
  """
  tt=model.get_travel_times(source_depth_in_km=depkm,distance_in_degree=distdeg,phase_list="P")[0].time
  return tt

def georpz(edist,edepth,p,alphar,betar,rhor,rearth):
  """
  Calculates the geometric spreading and surface response
  geometric spreading from Kanamori and Stewart (1976)
  surface excitation from Helmberger (1974)
  edist [degrees] and edepth [km] are for the event 
  alphar,betar,rhor are receiver values [km/s, gm/cm3]
  rearth in km
  returns g,rpz
  """
  #surface params
  salphar=6.5 #km/s
  srhor=2.9 #g/cm^3
  # sbetar=salphar/sqrt(3)  unused
  p=p*180/pi/rearth
  error=rearth/(rearth-edepth)
  angih=arcsin(p*salphar)
  #degih=angih*180/pi
  dd0=gttP(edepth,edist)-gttP(edepth,edist-1)
  dd1=gttP(edepth,edist+1)-gttP(edepth,edist)
  d2tddel2=dd1-dd0 # 2nd-order change in arrival time with distance (deg)
  pp=d2tddel2*(180/pi)**2
  dihdel=abs(pp/cos(angih)*salphar/(rearth-edepth))
  angi0=arcsin(p*alphar/error)
  factor=p*srhor*salphar**2/(rhor*alphar)
  g=(factor*dihdel/(sin(edist*pi/180)*cos(angi0)))**0.5
  #print("in georpz: ", factor,dihdel,edist,angi0,sin(edist*pi/180),cos(angi0))
  p2=(p/error)**2
  etaa=(1/(alphar**2)-p2)**.5
  etab=(1/(betar**2)-p2)**.5
  rpz=2*etaa*alphar*((etab**2-p2)/(betar**2*(etab**2-p2)**2+4*p2*betar**2*etaa*etab))
  #print("rpz params:",rpz,etaa,etab,p2,alphar,betar)
  return g,rpz

def estFgPcorrect(edistdeg):
  """
  Use Newman et al (1998) distance-based correction for dip-slip
  earthquakes. 
  Returns estFgP2
  """
  # polynomial spline fit to a number of real dip-slip EQs
  a0=1.17 ; a1=-7.27e-3 ; a2=6.01e-5
  estFgP2=a0 + a1*edistdeg + a2*edistdeg**2
  return estFgP2

def getFgP2(tr, Defaults=Defaults, Event=Event, FgP2min=0.2, **kwargs):
    trinfo=trstat2pd(tr)
    edistdeg=trinfo.iloc[0].distance[1]
    eqaz=trinfo.iloc[0].az
    edepth=Event.origin[0][2]
    phi,delta,lmbda=Event.focmech
    rho_site=Defaults.siteparams[1]
    pvel_site=Defaults.siteparams[0]
    rearth=Defaults.rearth
    aob=Defaults.aob
    qbc=Defaults.qbc
    p,g,rpz,fp,fpp,fsp,PP,SP=bc10(edistdeg,edepth,phi,delta,lmbda,eqaz,rho_site,pvel_site,rearth)
    FgP2calc=fp**2+(fpp*PP)**2+(2/(3*aob))*qbc*((SP*fsp)**2)
    if (FgP2calc<FgP2min):
        FgP2calc=FgP2min  # avoid blow-up
    # geometric spreading and near-surface excitation
    geomsp=g*rpz
    return FgP2calc,g,rpz,geomsp

def bc10(edist,edepth,phi,delta,lmbda,eqaz,rho,pvel,rearth):
  """
  Calculates geometric spreading, surface excitation, 
  radiation and reflection coefficients
  edist [degrees] and edepth [km] are for the event 
  phi, delta,lmbda is focal orientation (strike, dip,rake)
  pvel,rho are receiver values [km/s, gm/cm3]
  rearth in m (converted internally to km)
  
  from Boatwright & Choy (eqn 10) p.2097 
  returns p,g,rpz,fp,fpp,fsp,PP,SP
  """
  #edist  # in degrees
  #edepth # in km
  # near surface crustal params (input)
  # phif is phi relative to event (phi-takeoff)
  rhor=rho/1e3
  alphar=pvel/1e3
  betar=alphar/sqrt(3)
  rearth=rearth/1e3  #in km
  
  # rad coefficients
  d0=gttP(edepth,edist-0.5)
  d1=gttP(edepth,edist+0.5)
  dtddel=d1-d0  # change in arrival time with distance (deg)
  p=dtddel
  g,rpz = georpz(edist,edepth,p,alphar,betar,rhor,rearth)
  phif=phi-eqaz  # strike relative to earthquake azimuth
  # compute radiation coefficients
  aih=arcsin(p*alphar*180./(pi*rearth))
  #print("aih,alphar,p,pi,rearth",aih,alphar,p,pi,rearth)
  ajh=arcsin(sin(aih)/sqrt(3.))
  si=sin(aih)
  ci=cos(aih)
  #sj=sin(ajh)
  cj=cos(ajh)
  s2i=2.*si*ci
  #c2i=2.*ci*ci-1.
  sd=sin(delta*pi/180.)
  cd=cos(delta*pi/180.)
  sl=sin(lmbda*pi/180.)
  cl=cos(lmbda*pi/180.)
  sf=sin(phif*pi/180.)
  cf=cos(phif*pi/180.)
  s2f=2.*sf*cf
  c2f=2.*cf*cf-1.
  c2d=2.*cd*cd-1.
  sr=sd*cd*sl
  pr=cl*sd*s2f-sr*c2f
  qr=sl*c2d*sf+cl*cd*cf
  #print(delta,lmbda,phif,sd,cd,sl,cl,sf,cf)
  #pl=sr*s2f+cl*sd*c2f
  #ql=-cl*cd*sf+sl*c2d*cf
  fp=sr*(3.*ci*ci-1.)-qr*s2i-pr*si*si
  fpp=sr*(3.*ci*ci-1.)+qr*s2i-pr*si*si
  s2j=sin(2.*pi-2.*ajh)
  c2j=cos(2.*pi-2.*ajh)
  fsp=1.5*sr*s2j+qr*c2j+0.5*pr*s2j

  # reflection coefficients
  pbeta=p*betar*180./(pi*rearth)
  a=4.*pbeta**2*(betar/alphar)*ci*cj
  b=(1.-2.*pbeta*pbeta)**2.
  c=4.*(betar/alphar)*pbeta*cj*(1.-2.*pbeta**2)
  refPP=(a-b)/(a+b)
  refSP=c/(a+b)*ci/cj
  
  return p,g,rpz,fp,fpp,fsp,refPP,refSP

def gmeanCut(x,cutoff=0, **kwargs):
    """
    find geometric mean of the input array x
    If cutoff is supplied, it will remove values below or above *cutoff and /cutoff,
    otherwise it will give the straight geometric mean.  Program will also output new
    set of x-values  
    """
    x=x[np.nonzero(x)]  # only retain non-zero values (gmean will crash otherwise)
    xmean1=gmean(x)
    if cutoff>0:
        xkeep=x[x>xmean1/cutoff]
        xkeep=xkeep[xkeep<xmean1*cutoff]
        xmean=gmean(xkeep)
    else:
        xmean=xmean1
        xkeep=x
    return xmean,xkeep

def e2Me(e):
    """
    Outputs Energy Magnitude (Me) using the formalism of Choy and Boatright
    """
    return 2/3*log10(e)-2.9

class iterate:
    count = 0
    def str(self):
       return str(iterate.count).zfill(2)
    def step(self):
       iterate.count += 1


def eventdir(Defaults=Defaults,Event=Event,create=False,cd=True,**kwargs):
    """
    either go into an existing event directory, or create it and 
    then go in.
    Depending on what you're doing, it may be good to know where you are beforehand.
    """
    import os,shutil
    owd=os.getcwd()
    eyear=Event.origin[1].year
    edir=os.path.join(Defaults.edirbase,str(eyear),Event.eventname)  # event directory
    edirit=os.path.join(edir,Event.iter) # iteration path
    
    if create:   # will attempt to create a new directory        
        if os.path.exists(edirit):
            try:
                bakdir=edirit+".bak"
                shutil.move(edirit,bakdir)
            except:
                print("ERROR:  coudn't move directory to ",bakdir)
        try:
            os.makedirs(edirit)
            if cd:
                os.chdir(edirit)
                print("New working directory is: ",edirit)
        except:
            print("Error:  Couldn't create directory: ", edirit)
    else:   # use existing directory
        try:
            if cd:
                os.chdir(edirit)
                print("Working directory now: ",edirit)
        except:
            print("ERROR:  coudn't move into directory ", edirit)
    return edirit,owd

def getwaves(Defaults=Defaults, Event=Event, **kwargs): 
    """
    gets waveforms from NEIC for event using IRIS-supplied station information
    and store raw data in appropiate directory.

    outputs the event stream (st) and metadata datframe (df)
    """
    from obspy import UTCDateTime
    from obspy.core.stream import Stream
    from rtergpy.run import etime2name
    import os
    import pandas as pd
    import pickle

    eloc=Event.origin[0]
    etime=Event.origin[1]
    ecount=Event.ecount
    runiter=Event.iter
    eventname=Event.eventname
    network=Defaults.network
    chan=Defaults.chan
    rads=Defaults.stationrange
    pwindow=Defaults.waveparams[1]
    edirbase=Defaults.edirbase
    src=Defaults.src

    edirit,origwd=eventdir(Defaults=Defaults,Event=Event,create=True,cd=True)
    print("Checking for stations available within range from IRIS")
    inventory = get_respinv(network,eloc,etime,rads,chan) # from fsdn
    #print("Pulling Waveforms from NEIC")
    st = downloadwaves(inventory, eloc, etime, pwindow,src=src)  # stream
    if len(st) == 0:
        raise ValueError("ERROR: No waveforms obtained.") 
    now=UTCDateTime()
    # metadata I want to save for later. note, anything that is a list neds to be within brackets  
    df=pd.DataFrame({"eventname":eventname,"iteration":runiter,"etime":etime,"eloc":[eloc],
        "network":network,"chan":chan,"stationrange":[rads],"pwindow":[pwindow],
        "modtime":now,"eventdir":edirit,"inventory":[inventory]
        }, dtype=object)

    # create pkl directory
    ediritpkl=os.path.join(edirit, "pkls")
    if not os.path.exists(ediritpkl):
        os.mkdir(ediritpkl)
    # compression is fast and about 2x smaller
    # metadata
    dfpathfilepk=os.path.join(ediritpkl, "Params_"+eventname+".pkl")
    dfpathfilecsv=os.path.join(edirit, "Params_"+eventname+".csv")
    try:
        print("writing ",dfpathfilepk,"\n",dfpathfilecsv)
        df.to_pickle(dfpathfilepk)
        df.to_csv(dfpathfilecsv)
    except:
        print("ERROR:   writing",dfpathfilepk,"\n",dfpathfilecsv)

    # data raw and processesd
    stpathfile=os.path.join(ediritpkl, "Wavestream-raw_"+eventname+".pkl")
    try: 
        print("writing ",stpathfile)
        with open(stpathfile, 'wb') as file:
            pickle.dump(st,file)
    except:
        print("ERROR:   writing ",stpathfile)
    return st,df 

def loadwaves(Defaults=Defaults, Event=Event, **kwargs): 
    """
    loads existing waveform info
    """
    import pandas as pd
    eventname=Event.eventname
    edirbase=Defaults.edirbase
    runiter=str(Event.iter)
    eyr=Event.origin[1].year #eventname[0:4]   # get year
    edir=os.path.join(edirbase,str(eyr),eventname) # event directory
    ediritpkl=os.path.join(edir,runiter,"pkls")
    # load data 
    for file in os.listdir(ediritpkl):
        fpath=os.path.join(ediritpkl,file)
        # raw data load
        if file.startswith('Wavestream-raw_'+eventname):    
            if file.endswith('pkl'):
                st=pd.read_pickle(fpath)     
            elif file.endswith('pkl.gz'):
                st=pd.read_pickle(fpath, compression='gzip')     
        # load metadata
        if file.startswith("Params_"+eventname):
            df=pd.read_pickle(fpath)
    return st,df

def wave2energytinc(tr,Defaults=Defaults, Event=Event, fband=Defaults.waveparams[0][0], **kwargs):
  """
  calculate estimated earthquake energy as a function of time using 
  the rebuit time-series hopefully allowing for substantial speedup.
  frequency range, and earth parameters
  from Newman et al., (1998)
  tr = obspy style waveform with response information and distance info attached
  waveparams = [[fmin,fmax],[tstart,tlength,tstep=waveparams]]
  eqparams = [[elat,elon,edepth],[phi,delta,lmbda]]
  earthparams = [sitepvel,siterho]
  """
  from itertools import accumulate
  from matplotlib import pyplot as plt
  waveparams=Defaults.waveparams
  siteparams=Defaults.siteparams
  resample=Defaults.resample
  step=waveparams[2]
  
  fmin,fmax=fband
  prePtime,postPtime=waveparams[1]
  elat,elon,edepth=Event.origin[0]
  phi,delta,lmbda=Event.focmech
  pvel_site,rho_site=siteparams

  rearth=Defaults.rearth
  qbc=Defaults.qbc
  avfpsq=Defaults.avfpsq
  aob=Defaults.avfpsq 

  edistdeg=tr.stats.distdeg
  eqaz=tr.stats.az
  tr.resample(resample)
  t1cut=tr.stats.phasePtime+prePtime
  t2cut=tr.stats.phasePtime+postPtime
  trslice=tr.slice(t1cut,t2cut)
  tr=process_waves(trslice,freqmin=fmin,freqmax=fmax)
  trf=fft(tr)
  # # recreating obspy.freqatributes.spectrum
  n1=0; n2=len(tr)
  n=n2-n1
  srate=tr.stats.sampling_rate
  dt=1/srate    
  f=np.linspace(0,srate/2,n)

  # focal and distance corrections, 
  # focal corrections
  FgP2,g,rpz,geomsp=getFgP2(tr, Defaults=Defaults, Event=Event)
  # distance-based FgP2estimate
  estFgP2=estFgPcorrect(edistdeg)
  
  # integration prep 
  trftstar=np.full_like(trf.real,0)
  #sinu=0
  for j in range(0,len(f)-1):
    if (f[j]>fmin) & (f[j]<=fmax):
        trftstar[j]=trf[j].real*(exp(2*pi*f[j]*tstar(f[j])))**.5
    else:
        trftstar[j] = trf[j].real
  # remerge phase info with corrected amplitude and rebuild timeseries 
  ttstar=ifft(trftstar.real+1j*trf.imag) 
  correction1=2*pi*dt*rho_site*pvel_site*((rearth/geomsp)**2)*4*pi*avfpsq*(1+qbc) # true energy for mech
  #Ettstar=sum(abs(ttstar)**2)
  #Corrected_Energy=Ettstar*correction1/FgP2  # using supplied mechanism
  #Estimated_Energy=Ettstar*correction1/estFgP2  # averaged dip-slip
  # iterate over every value
  # energy per dt
  Ettstar_dt=list(accumulate([abs(i)**2 for i in ttstar]))
  Estimated_Energy_dt=[correction1/estFgP2*Edt for Edt in Ettstar_dt]
  est2corr=estFgP2/FgP2

  #report out by step (usually 1 second):
  stepn=int(srate*step)
  Estimated_Energy_tinc=(Estimated_Energy_dt[(stepn-1):])[::stepn]
  #print("Estimated/Corrected Energy_ps (last) = ", Estimated_Energy_tinc[-1],Corrected_Energy_tinc[-1])
  #plt.plot([i*dt for i in range(0,len(ttstar))],Estimated_Energy_tinc,'r-')
  return Estimated_Energy_tinc,estFgP2,FgP2,est2corr

def ErgsFromWaves(st,Defaults=Defaults,Event=Event,**kwargs):
    """
    Iterate through waves,time steps, 2 frequency bands to calculate energies
    returns 2 panda arrays of cumulative energy time series
    """
    waveparams=Defaults.waveparams
    siteparams=Defaults.siteparams
    fbands=waveparams[0]
    #tstart,tlength=waveparams[1]
    tstep=waveparams[2]
    elat,elon,edepth=Event.origin[0]
    phi,delta,lmbda=Event.focmech
    #pvel_site,rho_site=siteparams
    
    nsamples=waveparams[1][1]-waveparams[1][0]

    if len(fbands) > 2:
        print("WARNING:  will only iterate over first 2 of ",len(fbands),"frequency bands")
    fbandlabel="BB"
    for fband in fbands[:2]:   # will only iterate over first to sets of bands if more are included
        tempEdf=pd.DataFrame()  # energies
        print("Running fband",fband,"Hz:")
        netstatchan=[0]*len(st)
        fbandlist=[0]*len(st)
        waveparamlist=[0]*len(st)
        estFgP2=[0]*len(st)
        FgP2=[0]*len(st)
        est2corr=[0]*len(st)
        i = 0
        for tr in tqdm(st):
            # calc cum. energy and save in dataframe
            netstatchan[i]=str(tr).split(" | ")[0]
            Ergs=wave2energytinc(tr, Defaults, Event, fband=fband)
            # pad energy results with zeros for any waveforms that run short
            Epersec=list(Ergs[0])
            Epersec=(Epersec+nsamples *[0])[:nsamples]
            #print(i, netstatchan[i], len(Epersec))
            tempEdf[netstatchan[i]]=Epersec
            # build data frame with metadata, focal corrections
            fbandlist[i]=fband
            waveparamlist[i]=[waveparams[1],tstep]
            estFgP2[i]=Ergs[1]
            FgP2[i]=Ergs[2]
            est2corr[i]=Ergs[3]
            i += 1
        dfdict={"netstatchan":netstatchan,"fband"+fbandlabel:fbandlist,"waveparams":waveparamlist,
            "estFgP2":estFgP2,"FgP2":FgP2,"est2corr":est2corr}
        tempMDdf=pd.DataFrame(dfdict, dtype=object)
        if fbandlabel == "BB":   # BB first
            EBB=tempEdf
            Emd=tempMDdf
            fbandlabel="HF"
        else:     # then HF 
            EHF=tempEdf
            EHFmd=tempMDdf
    #Emd.rename(columns = {'fband':'fbandBB'}, inplace = True)
    Emd.insert(2, "fbandHF", EHFmd.fbandHF, True)
    return EBB,EHF,Emd
               #print(stationname,step,fband,thisE)

def trstat2pd(tr):
    """
    creates a pandas data frame from useful information within the trace.stats (and maybe other) sections
    retruns a data frame.
    """
    tr.stats
    dfdict={"netstatchan":str(tr).split(" | ")[0],
            "network":tr.stats.network,
            "station":tr.stats.station,
            "location":tr.stats.location,
            "channel":tr.stats.channel,
            "starttime":tr.stats.starttime,
            "endtime":tr.stats.endtime,
            "phasePtime": tr.stats.phasePtime,
            "sampling_rate":tr.stats.sampling_rate,
            "delta":tr.stats.delta,
            "npts":tr.stats.npts,
            "calib":tr.stats.calib,
            "az":tr.stats.az,
            "baz":tr.stats.baz,
            "coordinates": [[tr.stats.coordinates.latitude,tr.stats.coordinates.longitude]],
            "distance": [[tr.stats.distance,tr.stats.distdeg]],
            "pinc":tr.stats.pinc,
            "prayp":tr.stats.prayp,
            "ptoa":tr.stats.ptoa,
            }
    df=pd.DataFrame(dfdict, dtype=object)
    return df

def tacer(dE,prePtime=60, **kwargs):
    """
    Calculate the TACER from Convers and Newman 2013
    input is the time derivative of the cumulative energy  (dE)
    output is tacer (Time-Averaged Cumulative Energy Rate)
    """
    dEPonly=dE[int(0-prePtime):]
    tacerout=pd.DataFrame()
    itr=0  # trace number
    progressbar=tqdm(total=len(dEPonly.columns))
    while itr < len(dEPonly.columns):
        i=0  # position within trace
        cumerate=[]
        while i < len(dEPonly):
            i += 1
            cumerate.append(dEPonly.iloc[0:i,itr].sum()/(i))  # sum all values before ith locations
        tacerout[dEPonly.columns[itr]]=cumerate
        itr += 1
        progressbar.update(1)
    progressbar.close()
    return tacerout

def tacerstats(tacer):
    maxtacertime=tacer.idxmax()
    maxtacer=tacer.max()
    med=maxtacertime.median()
    m25=maxtacertime.quantile(0.25)
    m75=maxtacertime.quantile(0.75)
    # build df for return
    maxtacertime.name="time at max"
    maxtacer.name="max val"
    df=pd.merge(maxtacer,maxtacertime,right_index=True,left_index=True)

    return [med,m25,m75], df

### old  ###########################
 
def wave2energy(tr,waveparams,eqparams,siteparams):
  """
  Depricated? Using wave2energtinc now.
  calculate estimated earthquake for a waveform given earthquake parameters,
  frequency range, and earth parameters
  from Newman et al., (1998)
  tr = obspy style waveform with response information and distance info attached
  waveparams = [[fmin,fmax],[tstart,tlength,tstep=waveparams]]
  eqparams = [[elat,elon,edepth],[phi,delta,lmbda]]
  earthparams = [sitepvel,siterho]
  """
  fmin,fmax=waveparams[0]
  tstart,tlength=waveparams[1]
  elat,elon,edepth=eqparams[0]
  phi,delta,lmbda=eqparams[1]
  pvel_site,rho_site=siteparams

  rearth=6371e3  # meters
  qbc=15.6 #q-factor from B&C
  avfpsq=(4/15)
  aob=3**.5  # alpha over beta

  edistdeg=tr.stats.distdeg
  eqaz=tr.stats.az
  tr.resample(resample)
  t1cut=tr.stats.phasePtime-tstart 
  t2cut=t1cut+tlength
  trslice=tr.slice(t1cut,t2cut)
  tr=process_waves(trslice,freqmin=fmin,freqmax=fmax)
  #tr=trslice
  trf=fft(tr)
  # # recreating obspy.freqatributes.spectrum
  n1=0; n2=len(tr)
  n=n2-n1
  srate=tr.stats.sampling_rate
  dt=1/srate    
  f=np.linspace(0,srate/2,n)
  
  # for testing run the following in the data directory:
  # printf "II.WRAB_.00.BHZ.SAC\nII.WRAB_.00.BHZ.SAC.pz\n300\n1\n"  | /home/jconvers/EQerg/CWBerg/src/nergy_wprep.cwb.wderivs_outputfft
  # integration prep 
  sinu=0
  for j in range(0,len(f)-1):
    if (f[j]>fmin) & (f[j]<=fmax):
      sinu += abs(trf[j])**2*exp(2*pi*f[j]*tstar(f[j]))
  sinu=sinu*2*pi*dt/n   # not 100% on the 2pi
  Estar=rho_site*pvel_site*sinu  # didn't divide by pi, as this isn't called for in discrete Parservel's theorem
  p,g,rpz,fp,fpp,fsp,PP,SP=bc10(edistdeg,edepth,phi,delta,lmbda,eqaz,rho_site,pvel_site,rearth)
  geomsp=g*rpz
  #print("g=",g)
  #print("rpz",rpz)
  Energy=Estar*(rearth/geomsp)**2
  FgP2=fp**2+(fpp*PP)**2+(2/(3*aob))*qbc*((SP*fsp)**2)
  if (FgP2<0.2):
    FgP2=0.2  # avoid blow-up
  Nergy=Energy*4*pi*(avfpsq/FgP2)*(1+qbc) # true energy for mech
  # distance-based FgP2estimate
  estFgP2=estFgPcorrect(edistdeg)
  estNergy=Energy*4*pi*(avfpsq/estFgP2)*(1+qbc) # true energy for mech
  #print("sinu,Estar,Energy = ", sinu,Estar,Energy,Nergy,estNergy)
  return estNergy,Nergy,estFgP2,FgP2
