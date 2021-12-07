from obspy.core.stream import Stream
from obspy import UTCDateTime
from numpy import sin,cos,arcsin,sqrt,abs,pi,log10,exp
import os
import numpy as np
import pandas as pd
from rtergpy.waveforms import gmeanCut,tacerstats
from rtergpy.run import defaults
Defaults=defaults()
import matplotlib.pyplot as plt

def fbandlabels(Emd):
    fbb=Emd.fbandBB[0]
    fhf=Emd.fbandHF[0]
    if  fbb[0] < 0.00001:
        fminbblabel="0"
    else :
        fminbblabel=str("1/"+str(int(1/fbb[0])))
    labelbb=fminbblabel+" - "+str(int(fbb[1]))+" Hz"
    labelhf="1/"+str(int(1/fhf[0]))+" - "+str(int(fhf[1]))+" Hz"
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
    while i < len(tacer.columns):
        az=float(trdf[trdf['netstatchan'] == tacer.columns[i]]['az']) # azimuth (need to search since some are dropped) 
        [dist]=list(trdf['distance'][trdf['netstatchan'] == tacer.columns[i]])
        plt.plot(amp*tacer.iloc[:,i]/max(tacer.iloc[:,i])+az-amp, '-', color='red', lw=0.8, alpha=0.6) 
        plt.plot(meds.iloc[i,1],az, 'r|', alpha=1) 
        plt.plot(meds.iloc[i,1],az, 'r.', alpha=1) 
        i +=1 
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

    labelbb,labelhf=fbandlabels(Emd)
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    plt.plot(distance,ebb/emeanbb,'ro',alpha=0.6,label=labelbb)
    plt.plot(distance,ehf/emeanhf,'bo',alpha=0.6,label=labelhf)
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

    labelbb,labelhf=fbandlabels(Emd)
    fig=plt.figure()
    fig.patch.set_facecolor('white')
    plt.plot(trinfo["az"],ebb/emeanbb,'ro',alpha=0.6,label=labelbb)
    plt.plot(trinfo["az"],ehf/emeanhf,'bo',alpha=0.6,label=labelhf)
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

    labelbb,labelhf=fbandlabels(Emd)
    #print("%s: %d traces, %.2e +- 10^%.2f [J]" %(labelbb,len(keepbb), emeanbb, np.std(np.log10(keepbb))))
    #print("%s: %d traces, %.2e +- 10^%.2f [J]" %(labelhf,len(keephf), emeanhf, np.std(np.log10(keephf))))

    fig=plt.figure()
    fig.patch.set_facecolor('white')
    # histogram it 
    #kwargs=dict(histtype="stepfilled", alpha=0.6, bins=10)
    kwargs=dict(histtype="bar", alpha=0.6, bins=9)   
    plt.hist(np.log10(keepbb/emeanbb), color='red', label=labelbb, **kwargs)
    plt.hist(np.log10(keephf/emeanhf), color='blue', label=labelhf, **kwargs)
    plt.legend()
    plt.ylabel('Count')
    plt.xlabel('Per station energy relative to mean (log_10)')
    plt.title(eventname+' Spread in energy results per station')
    plt.savefig('Ehist_'+eventname+'.png', dpi=150)
    #print(labelhf, emeanhf,len(keephf), np.std(np.log10(keephf)))
    if show:
        plt.show()
 
def Etincplot(Ebb,Ehf,Emd,trdf,ntr=0, **kwargs):
    labelbb,labelhf=fbandlabels(Emd)
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
    mapcpt=os.path.join(libdir,'map_gray.cpt')
    with gmt.config(MAP_FRAME_PEN="0.3p,black"):
        fig.grdimage(grid=grid, cmap=mapcpt, shading=True, projection=proj)
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
    enorm=np.log10(eint/emean)

    lats=[] ; lons=[]
    for netstatchan in E.columns:
        coords= trdf[trdf['netstatchan'] == netstatchan]['coordinates'].tolist()[0]
        lons.append(coords[1])
        lats.append(coords[0])

# plot stations with relative energy
    gmt.makecpt(cmap="split",background="True", continuous="True", series=[-1,1,0.5])
    fig.plot(x=lons,y=lats,style='t0.5c', color=enorm, cmap=True, pen="0.5,white",transparency=20)
    fig.colorbar(position="13.7/0.2+w2.5c/0.2c+v")

# plot event location in center
    fig.plot(x=eloc[1],y=eloc[0], style="a0.85c",color="green",pen="1p,black", transparency=50)

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