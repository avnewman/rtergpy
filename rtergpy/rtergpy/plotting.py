from obspy.core.stream import Stream
from obspy import UTCDateTime
from numpy import sin,cos,arcsin,sqrt,abs,pi,log10,exp
import os
import numpy as np
import pandas as pd
from rtergpy.waveforms import gmeanCut,e2Me #,tacerstats
from rtergpy.run import defaults,event
Event=event()
Defaults=defaults()
import matplotlib.pyplot as plt

def fbandlabels(wp0):
    fbb=wp0[0]
    fhf=wp0[1]
    
    if  fbb[0] < 0.00001:
        fminbblabel="0"
    else :
        fminbblabel=str("1/"+str(int(1/fbb[0])))
    if fbb[1] < 1 :
        fmaxbblabel=str("1/"+str(int(1/fbb[1])))
    else:
        fmaxbblabel=str(int(fbb[1]))
    
    if  fhf[0] < 0.00001:
        fminhflabel="0"
    elif fhf[0] < 1 :
        fminhflabel=str("1/"+str(int(1/fhf[0])))
    else :
        fminhflabel=str(int(fhf[0]))
    if fhf[1] < 1 :
        fmaxhflabel=str("1/"+str(int(1/fhf[1])))
    else:
        fmaxhflabel=str(int(fhf[1]))
    
    labelbb=fminbblabel+" - "+fmaxbblabel+" Hz"
    labelhf=fminhflabel+" - "+fmaxhflabel+" Hz"
    return labelbb,labelhf

def tacerplot(tacer,trdf,ttimes,meds,eventname,amp=50,show=True, **kwargs):
    """
    plots tacer as a function of azimuth, removing values below the SNR ratio threshold, but showing all.
    Will output the median, 25 and 75 percentile results.
    """

    i=0
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    plt.axvline(x=ttimes[0],color='red', ls='-',alpha=0.6,lw=2)
    npts=len(tacer.iloc[:,0])
    normsum=pd.Series(np.zeros(npts))
    
    while i < len(tacer.columns):
        az=float(trdf[trdf['netstatchan'] == tacer.columns[i]]['az']) # azimuth (need to search since some are dropped) 
        [dist]=list(trdf['distance'][trdf['netstatchan'] == tacer.columns[i]])
        norm=tacer.iloc[:,i]/max(tacer.iloc[:,i])
        plt.plot(amp*(norm-1)+az, '-', color='red', lw=0.8, alpha=0.6) 
        plt.plot(meds.iloc[i,1],az, 'r|', alpha=1) 
        plt.plot(meds.iloc[i,1],az, 'r.', alpha=1) 
        normsum += norm
        i +=1 
    normsumnorm=normsum/max(normsum)
    (normsumnorm*360-amp).plot(color='black', lw=1, alpha=0.5)
    
    plt.ylim(-amp,360)
    plt.xlim(-5,min([ttimes[2]*3,len(tacer)]))
    plt.xlabel("Time relative to predicted P-arrival [s]")
    plt.ylabel("Station azimuth from event (at peak)")
    plt.title(eventname+" TACER with azimuth")
    timing=str(ttimes[0])+" -/+ "+str(ttimes[1])+"/"+str(ttimes[2])+" s (25/75th pcnt.)"
    plt.text(0.88,0.84,timing, fontsize=9, horizontalalignment='right',verticalalignment='center', transform=fig.transFigure )
    plt.axvline(x=0,color='black', ls='--',lw=1)
    plt.savefig('Tacer_'+eventname+'.png', dpi=150)
    if show:
        plt.show()
 
def Edistplot(Ebb,Ehf,Emd,trinfo,eventname,ttime, prePtime=-60,cutoff=15,show=True, **kwargs):
    intval=int(ttime-prePtime)
    ebb=np.array(Ebb.iloc[intval])
    ehf=np.array(Ehf.iloc[intval])
    distance=np.array(trinfo["distance"].tolist())[:,1]   # in degree (0 is in meters)
    emeanbb,keepbb=gmeanCut(ebb,cutoff=cutoff)
    emeanhf,keephf=gmeanCut(ehf,cutoff=cutoff) 

    labelbb,labelhf=fbandlabels([Emd.fbandBB[0],Emd.fbandHF[0]])
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    plt.plot(distance,ebb/emeanbb,'o',alpha=0.6,label=labelbb, color='blue')
    plt.plot(distance,ehf/emeanhf,'o',alpha=0.6,label=labelhf, color='red')
    plt.plot(distance,[1 for i in distance],'k-') # line at 1
    plt.yscale('log')
    plt.legend()
    plt.xlabel('Distance [degrees]')
    plt.ylabel('Normalized Energy at '+str(ttime)+'s [J]')
    plt.title(eventname+' Station energy with distance')
    plt.ylim(.1,10)
    plt.savefig('Edist_'+eventname+'.png', dpi=150)
    if show:
        plt.show()
 
def Eazplot(Ebb,Ehf,Emd,trinfo,eventname,ttime, prePtime=-60,cutoff=15,show=True, **kwargs):
    intval=int(ttime-prePtime)
    ebb=np.array(Ebb.iloc[intval])
    ehf=np.array(Ehf.iloc[intval])
    emeanbb,keepbb=gmeanCut(ebb,cutoff=cutoff)
    emeanhf,keephf=gmeanCut(ehf,cutoff=cutoff)

    labelbb,labelhf=fbandlabels([Emd.fbandBB[0],Emd.fbandHF[0]])
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    plt.plot(trinfo["az"],ebb/emeanbb,'o',alpha=0.6,label=labelbb, color='blue')
    plt.plot(trinfo["az"],ehf/emeanhf,'o',alpha=0.6,label=labelhf, color='red')
    plt.plot(trinfo["az"],[1 for i in trinfo["az"]],'k-')
    plt.legend()
    plt.xlabel('Azimuth [degrees]')
    plt.ylabel('Normalized Energy at '+str(ttime)+'s [J]')
    plt.title(eventname+' Station energy with azimuth')
    plt.yscale('log')
    plt.ylim(.1,10)
    plt.savefig('Eaz_'+eventname+'.png',dpi=150)
    if show:
        plt.show()
    
def Ehistogram(Ebb,Ehf,Emd,eventname,ttime, prePtime=-60,cutoff=15,show=True, **kwargs):
    intval=int(ttime-prePtime)
    ebb=np.array(Ebb.iloc[intval])
    ehf=np.array(Ehf.iloc[intval])
    emeanbb,keepbb=gmeanCut(ebb,cutoff=cutoff)
    emeanhf,keephf=gmeanCut(ehf,cutoff=cutoff)

    labelbb,labelhf=fbandlabels([Emd.fbandBB[0],Emd.fbandHF[0]])

    fig=plt.figure()
    fig.patch.set_facecolor('white')
    # histogram it 
    #kwargs=dict(histtype="stepfilled", alpha=0.6, bins=10)
    kwargs=dict(histtype="bar", alpha=0.6, bins=9)   
    plt.hist(np.log10(keepbb/emeanbb), color='blue', label=labelbb, **kwargs)
    plt.hist(np.log10(keephf/emeanhf), color='red', label=labelhf, **kwargs)
    plt.legend()
    plt.ylabel('Count')
    plt.xlabel('Per station energy relative to mean (log_10)')
    plt.title(eventname+' Spread in energy results per station')
    plt.savefig('Ehist_'+eventname+'.png', dpi=150)
    #print(labelhf, emeanhf,len(keephf), np.std(np.log10(keephf)))
    if show:
        plt.show()
 
def Etincplot(Ebb,Ehf,Emd,trdf,ntr=0, **kwargs):
    labelbb,labelhf=fbandlabels([Emd.fbandBB[0],Emd.fbandHF[0]])
    etincBB=np.array(Ebb.iloc[:,ntr])
    etincHF=np.array(Ehf.iloc[:,ntr])
    
    ax1=plt.gca()
    ax2=plt.twinx()

    ax1.plot(etincBB,'k-',alpha=0.6, label=labelbb)
    ax2.plot(etincHF,'b-',alpha=0.4, label=labelhf)
    ax1.set_ylabel(str("BB Energy Growth ["+labelbb+"]"))
    ax2.set_ylabel(str("HF Energy Growth ["+labelhf+"]"),color='blue')

def stationEmapPygmt(E,eloc,trdf,eventname,ttime,prePtime=-60,cutoff=15,itername=0,libdir=Defaults.libdir,show=True,**kwargs):
    import pygmt as gmt 
    import os
    intval=int(ttime-prePtime)   
    fig=gmt.Figure()

# make basemap    
    scale=15
    proj="A"+str(eloc[1])+"/"+str(eloc[0])+"/"+str(scale)+"c"
    grid = gmt.datasets.load_earth_relief(resolution="30m")
#    dgrid = gmt.grdgradient(grid=grid, radiance=[45, 80])
    mapcpt=os.path.join(libdir,'map_gray.cpt')
    with gmt.config(MAP_FRAME_PEN="0.3p,black"):
        fig.grdimage(grid=grid, cmap=mapcpt, shading=True, projection=proj)
#        fig.grdimage(grid=dgrid, cmap=mapcpt, projection=proj)
        fig.coast(shorelines="0.2p,50", area_thresh=100000, frame="afg30", projection=proj)

# background earthquakes
    eqcpt=os.path.join(libdir,'eq_color.cpt')
    eqs=os.path.join(libdir,'EQM4plus.data.sorted')
    fig.plot(data=eqs, style='c', cmap=eqcpt, pen='0.2,100', transparency=50)

# Plates
    plates=os.path.join(libdir,'PB2002_plates.txt')
    fig.plot(data=plates, pen='0.5,0', transparency=10)

# radial lines (approximate)
    degscale=220 # I don't know why this is the value, but works
# thick lines at distance limits
    with open('rads.txt','w') as file:
        file.write(str(eloc[1])+" "+str(eloc[0])+" "+str(25*degscale)+"\n")
        file.write(str(eloc[1])+" "+str(eloc[0])+" "+str(80*degscale)+"\n")
    fig.plot(data='rads.txt',style='E-',pen="1.5p,20",projection=proj)
 
# every 20 degrees dashed line
    with open('rads20.txt','w') as file:
        for dist in range(10,80,10):
            file.write(str(eloc[1])+" "+str(eloc[0])+" "+str(dist*degscale)+"\n")
    fig.plot(data='rads20.txt',style='E-',pen="0.5p,50,-",projection=proj)
    os.remove('rads.txt')
    os.remove('rads20.txt')
# prep data for plotting
    eint=np.array(E.iloc[intval])
    emean,keep=gmeanCut(eint,cutoff=cutoff)
    #print(eint)   # need to look at, giving issues with division by zero
    enorm=np.log10(eint/emean)

    lats=[] ; lons=[]
    for netstatchan in E.columns:
        coords= trdf[trdf['netstatchan'] == netstatchan]['coordinates'].tolist()[0]
        lons.append(coords[1])
        lats.append(coords[0])

# plot stations with relative energy
    #gmt.makecpt(cmap="split",background="True", continuous="True", series=[-1,1,0.5])
    gmt.makecpt(cmap="split",background="True", series=[-1,1,0.25])
    fig.plot(x=lons,y=lats,style='t0.5c', fill=enorm, cmap=True, pen="0.5,white",transparency=20)
    fig.colorbar(position="x13.7/0.2+w2.5c/0.2c+v")

# plot event location in center
    fig.plot(x=eloc[1],y=eloc[0], style="a0.85c",fill="green",pen="1p,black", transparency=50)

# add text
    with gmt.config(MAP_FRAME_PEN="0,white"):
        proj2="X"+str(scale)+"c"
        fig.basemap(region=[0,1,0,1], projection=proj2, frame="wesn")
        fig.text(text='Log@-10@- (@~D@~ Energy)',x=0.89,y=0.015, justify="LM", font="9p,Helvetica,black", angle=90)
        fig.text(text='Iteration: '+str(itername),x=0.87,y=0.97, justify="LM", font="10p,Helvetica,black", angle=0)
        fig.text(text='Event: '+str(eventname),x=0.015,y=0.97, justify="LM", font="10p,Helvetica,black", angle=0)
        for deg in [25,80]:
            xpos=-(.707*deg/90/2)+0.5-0.006 # SW quad 45deg x=y
            fig.text(text=str(deg)+'@.',x=xpos,y=xpos, justify="CM",fill=220, font="10p,Helvetica,black", angle=-45)

    fig.savefig('StationMap_'+eventname+'.png', dpi=150)
    if show:
        fig.show()

def fracpostID(df,frac=0.5, **kwargs):
    iMax=df[df.eq(1)].index[0]
    PostiMax=df.iloc[iMax:] 
    ltfrac=PostiMax[PostiMax.lt(frac)]
    return ltfrac.index[0]

def Efluxplots(dEdt, trdf, Event=event(), Defaults=defaults(), ampDIST=5, ampAZ=30, pcts=[0.5,0.25,0.1], show=True, **kwargs):
    from obspy import UTCDateTime
    #UTCDateTime.DEFAULT_PRECISION=1

    eventname=Event.eventname
    eloc, etime=Event.origin
    prePtime,postPtime=Defaults.waveparams[1]
    nst=len(dEdt.iloc[0])  # N stations
    npts=postPtime-prePtime # N time-series points
    xmax=postPtime
 
 # distance plot (top)
    plt.subplots(figsize=(12,12), facecolor='white')
    plt.subplot(3, 1,1)
    
    normsum=pd.Series(np.zeros(npts))
    Deltamin,Deltamax=Defaults.stationrange

    for count in range(nst):
        delta=float(trdf.distance[count][1])
        max=np.max(dEdt.iloc[:,count])
        norm=(dEdt.iloc[:,count]/max).clip(lower=0) # don't allow negative values to creep in for missing/strange data 
        (norm.shift(prePtime,axis=0)*ampDIST+delta).plot(color='gray',alpha=0.5)
        normsum += norm.clip(lower=0)  
    normsummax=np.max(normsum)
    normsumnorm=normsum/normsummax
    #(Deltamin+normsumnorm.shift(prePtime,axis=0)*(Deltamax-Deltamin)).plot(color='black')
    plt.axvline(x=0,color='black', ls='--', lw=1)
    plt.ylim(Deltamin,Deltamax+ampDIST)
    plt.title('Smoothed Energy Flux (normalized and stacked) [J/s]')
    #plt.xlabel('Time from Theor. P-arrival [s]')
    plt.gca().xaxis.grid(True)
    plt.ylabel('Distance [°]')
   
    timestr="%4d/%02d/%02d %02d:%02d:%04.1f" % (etime.year,etime.month,etime.day,etime.hour,etime.minute,etime.second)
    locstr="%.2f°, %.2f°, %.1f km" % (eloc[0],eloc[1],eloc[2])
    magstr="Mw=%.1f, Me=%.1f, T_tac=%d s" % (Event.Mw, Event.Me, Event.ttime)
    plt.text(xmax,(Deltamax+ampDIST)*.98, magstr +': '+timestr+' @ '+locstr, 
        va='top', ha='right', size='medium',
        bbox=dict(boxstyle='round',facecolor='white', alpha=0.8))

    # azimuth plot
    plt.subplot(3, 1,2)
    for count in range(nst):
        az=float(trdf.az[count])
        max=np.max(dEdt.iloc[:,count])
        norm=(dEdt.iloc[:,count]/max).clip(lower=0) # avoid negs as above
        (norm.shift(prePtime,axis=0)*ampAZ+az).plot(color='gray',alpha=0.5)
    plt.axvline(x=0,color='black', ls='--', lw=1)
    #(normsumnorm.shift(prePtime,axis=0)*360).plot(color='black')
    plt.ylim(0,360+ampAZ)
 #   plt.xlabel('Time from Theor. P-arrival [s]')
    plt.ylabel('Azimuth [°]')
    plt.tight_layout()
    plt.gca().xaxis.grid(True)
    
    # stacked solution
    plt.subplot(3,1,3)
    for count in range(nst):
        az=float(trdf.az[count])
        max=np.max(dEdt.iloc[:,count])
        norm=dEdt.iloc[:,count]/max
        (norm.shift(prePtime,axis=0)).plot(color='gray',alpha=0.3)
        
    plt.axvline(x=0,color='black', ls='--', lw=1)
    (normsumnorm.shift(prePtime,axis=0)).plot(color='black')

    #pcts=[0.1,0.2,0.5] # list of fractions
    droptimes=[]
    for pct in pcts:
        droptime=fracpostID(normsumnorm,frac=pct)+prePtime
        droptimes.append(droptime)
        plt.axvline(x=droptime,color='black', ls='-.', lw=1 )
        plt.text(droptime,1, str(droptime)+' @ '+str(pct)+' max', rotation=90, va='top', ha='right', bbox=dict(boxstyle='round',facecolor='white', alpha=0.4))

    plt.ylim(0,1.05)
    plt.xlabel('Time from Theor. P-arrival [s]')
    plt.ylabel('Stacked E-flux')
    plt.tight_layout()
 
    plt.gca().xaxis.grid(True)
    plt.savefig('Eflux_Dist-Az_'+eventname+'.png', dpi=150)
    if show:
        plt.show()
    return droptimes

def Efluxplots_old(dEdt,trdf,eventname, prePtime=-60, ampDIST=10, ampAZ=30, show=True, **kwargs):
    npts=len(dEdt.iloc[:,0]) # N time-series points
    nst=len(dEdt.iloc[0])  # N stations
    normsum=pd.Series(np.zeros(npts))
    for count in range(nst):
        delta=float(trdf.distance[count][1])
        max=np.max(dEdt.iloc[:,count])
        norm=dEdt.iloc[:,count]/max
        (norm.shift(prePtime,axis=0)*ampDIST+delta).plot(color='gray',alpha=0.5)
        normsum += norm 
    normsummax=np.max(normsum)
    normsumnorm=normsum/normsummax
    (normsumnorm.shift(prePtime,axis=0)*80).plot(color='black')
    plt.axvline(x=0,color='black', ls='--', lw=1)
    plt.ylim(0,80+ampDIST)
    #plt.xlim(0,600)
    plt.title('Smoothed Energy Flux (normalized and stacked) [J/s]')
    plt.xlabel('Time from Theor. P-arrival [s]')
    plt.ylabel('Distance [°]')
    plt.savefig('Eflux_Dist_'+eventname+'.png', dpi=150)
    if show:
        plt.show()

    # azimuth plot
    for count in range(nst):
        az=float(trdf.az[count])
        max=np.max(dEdt.iloc[:,count])
        norm=dEdt.iloc[:,count]/max
        (norm.shift(prePtime,axis=0)*ampAZ+az).plot(color='gray',alpha=0.5)
        #print(count,delta,np.max(normsum))
    #plt.show()
    (normsumnorm.shift(prePtime,axis=0)*360).plot(color='black')
    plt.axvline(x=0,color='black', ls='--', lw=1)
    plt.ylim(0,360+ampAZ)
    #plt.xlim(60,600)
    plt.title('Smoothed Energy Flux (normalized and stacked) [J/s]')
    plt.xlabel('Time from Theor. P-arrival [s]')
    plt.ylabel('Azimuth [°]')
    plt.savefig('Eflux_Az_'+eventname+'.png', dpi=150)
    if show:
        plt.show()


### Below are testing 
def stationEmapBasemap(E,eloc,trdf,intval=-1,cutoff=15,**kwargs):
    import matplotlib as mpl
    from mpl_toolkits.basemap import Basemap

    eint=np.array(E.iloc[intval])
    emean,keep=gmeanCut(eint,cutoff=cutoff)
    lats=[] ; lons=[]
    for netstatchan in E.columns:
        coords= trdf[trdf['netstatchan'] == netstatchan]['coordinates'].tolist()[0]
        lons.append(coords[1])
        lats.append(coords[0])

    enorm=np.log10(eint/emean)
    
    plt.rcParams['figure.dpi'] = 250  # this controls the size
    plt.figure(figsize=(20,20))
    map = Basemap(projection='ortho',lat_0=eloc[0],lon_0=eloc[1], resolution='l')
    map.drawcoastlines(linewidth=0.25)
    map.drawcountries(linewidth=0.25)
    map.drawmapboundary(fill_color='whitesmoke')
    map.drawmeridians(np.arange(0,360,20))
    map.drawparallels(np.arange(-90,90,30))
    map.bluemarble(scale=.1,alpha=0.45)
    x,y=map(lons,lats)
    #ax.stock_img()  # background (wish I could make it B&W)
    cpt=plt.cm.seismic
   
    plt.scatter(x,y,c=enorm,cmap=cpt,norm=mpl.colors.Normalize(vmin=-1.5,vmax=1.5),marker='^', edgecolors='black', s=800, alpha=0.8)
    ex,ey=map(eloc[1],eloc[0])
    plt.plot(ex,ey, marker='*', color='y', ms=60, mec='k', alpha=0.6)

    plt.show()

def ETxoplot(EBBlmean,EHFlmean,upFitResults,downFitResults,Event=event(),Defaults=defaults(),show=False):

    # Extract best window fit info
    min1Window, min1Peak_t1, min1Peak_t2, min1Peak_slope, min1Peak_intercept, min1MisfitPeak, m1results = upFitResults
    min2Window, min2Peak_t1, min2Peak_t2, min2Peak_slope, min2Peak_intercept, min2MisfitPeak, m2results = downFitResults
    prePtime=Defaults.waveparams[1][0]
    EBBTxo=Event.EBBTxo
    EHFTxo=Event.EHFTxo
    Txo=Event.Txo
    labelbb,labelhf=fbandlabels(Defaults.waveparams[0])
    HFEcorr=Defaults.HFEcorr

    # Generate x and y for the fit plots
    x=np.arange(0, len(EHFlmean))
    y1_fitFull= min1Peak_slope * x + min1Peak_intercept
    y2_fitFull= min2Peak_slope * x + min2Peak_intercept
    x1_fit = np.arange(min1Peak_t1, min1Peak_t2)
    y1_fit = min1Peak_slope * x1_fit + min1Peak_intercept
    x2_fit = np.arange(min2Peak_t1, min2Peak_t2)
    y2_fit = min2Peak_slope * x2_fit + min2Peak_intercept

    # Plot the data and the fit
    tsize=len(EBBlmean)+prePtime
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, hspace=0.1, height_ratios=[2, 2])
    axs = gs.subplots(sharex=True, sharey=False)
    fig.suptitle(f"Cumulative Energy Growth ($M_ETxo$= {e2Me(EBBTxo):.2f}, Txo= {Txo:.1f} s)", fontsize=16, y=0.965)

    # top plot
    axs[0].plot(x+prePtime, EBBlmean, label="",color='blue', linewidth=3)
    axs[0].set_ylabel("Cumulative BB Energy")
    axs[0].axvline(x=Txo, color='black', linewidth=3, linestyle='-', alpha=1, label=f'Txo = {Txo:.2f}s')
    axs[0].set_ylim(0, 1.1* EBBlmean.max())
    axs[0].axhline(y=EBBTxo, color='black', linewidth=1, linestyle=':', alpha=0.5, label='')
    axs[0].text(tsize*0.95, 1.05*EBBTxo, f'$E_BB(Txo)$ = {EBBTxo:.1e} J, $M_E$ = {e2Me(EBBTxo):.2f}', horizontalalignment='right')
    axs[0].text(tsize*0.95,0.09*EBBlmean.max(), f"N stats used = {Event.nBBlmean}", horizontalalignment='right')
    axs[0].text(tsize*0.95,0.03*EBBlmean.max(), f"f-band ({labelbb})", horizontalalignment='right')
    
    title=Event.origin[1].strftime("%Y/%m/%d %H:%M:%S")+" ("+str(Event.eventname)+f") at {Event.origin[0][0]:.2f}°, {Event.origin[0][1]:.2f}°, Z={Event.origin[0][2]:.1f} "
    axs[0].set_title(title)

    # right side of plot
    ax0r = axs[0].twinx()
    MeGrowth = e2Me(EBBlmean)

    # create ticks for Me axis
    start = round(MeGrowth.max()-0.7, 1)
    end = round(MeGrowth.max(), 1)
    tick_labels = [round(val, 1) for val in np.arange(start, end + 0.1, 0.1)]
    Ebb_ticks = [10 ** ((me + 2.9) * 1.5) for me in tick_labels]
    ax0r.set_yticks(Ebb_ticks)
    ax0r.set_yticklabels(tick_labels)
    ax0r.set_ylabel("$M_Ebb$ (BB Energy Magnitude)")
    ax0r.tick_params(axis='y', colors='black')
    ax0r.spines['right'].set_color('black')
    ax0r.set_ylim(0, 1.1* EBBlmean.max())  # set limits for Me axis

    # bottom plot
    axs[1].plot(x+prePtime, EHFlmean, label="", color='red',linewidth=3)  # average of all tacer values
    axs[1].plot(x+prePtime, y1_fitFull, '-', label="",color='gray', linewidth=1,alpha=0.5)
    axs[1].plot(x+prePtime, y2_fitFull, '-', label="Best fits", color='gray', linewidth=1,alpha=0.5)
    axs[1].plot(x1_fit+prePtime, y1_fit, '-', label="",color='black', linewidth=3, alpha=0.5)
    axs[1].plot(x2_fit+prePtime, y2_fit, '-', label="", color='black', linewidth=3, alpha=0.5)
    axs[1].set_ylim(0, 1.1* EHFlmean.max())
    axs[1].set_xlim(0, tsize)
    axs[1].axhline(y=EHFTxo, color='black', linewidth=1, linestyle=':', alpha=0.5, label=f'ETxo = {EHFTxo:.1e}s')
    axs[1].text(tsize*0.95, 1.05*EHFTxo, f'$E(Txo)$ = {EHFTxo:.1e} J, $M_E$ = {e2Me(EHFTxo,eCorrection=HFEcorr):.2f}', horizontalalignment='right')

    axs[1].axvline(x=Txo, color='black', linewidth=3, linestyle='-', alpha=1, label=f'Txo = {Txo:.2f}s')
    axs[1].set_xlabel("Time from Theoretical P-arrival [s]")
    axs[1].set_ylabel("Cumulative HF Energy")
    axs[1].text(tsize*0.95,0.09*EHFlmean.max(), f"N stats used = {Event.nHFlmean}", horizontalalignment='right')
    axs[1].text(tsize*0.95,0.03*EHFlmean.max(), f"f-band ({labelhf})", horizontalalignment='right')
    
    axs[1].text(Txo+5,0.05*EHFlmean.max(), f"Txo= {Txo:.1f} s", horizontalalignment='left', rotation=90)

    # right side of plot
    ax1r = axs[1].twinx()
    MeGrowth = e2Me(EHFlmean, eCorrection=HFEcorr)

    # create ticks for Me axis
    start = round(MeGrowth.max()-0.7, 1)
    end = round(MeGrowth.max(), 1)
    tick_labels = [round(val, 1) for val in np.arange(start, end + 0.1, 0.1)]
    Ehf_ticks = [10 ** ((me + 2.9) * 1.5) for me in tick_labels]
    Ehf_ticks = np.array(Ehf_ticks)/HFEcorr  # correct for HFEcorr
    ax1r.set_yticks(Ehf_ticks)
    ax1r.set_yticklabels(tick_labels)
    ax1r.set_ylabel("$M_Ehf$ (HF Energy Magnitude)")
    ax1r.tick_params(axis='y', colors='black')
    ax1r.spines['right'].set_color('black')
    ax1r.set_ylim(0, 1.1* EHFlmean.max())  # set limits for Me axis
    plt.savefig('ETxo_'+Event.eventname+'.png', dpi=150)
    if show:
        plt.show()