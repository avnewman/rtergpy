#!/usr/bin/env python
# Script to calculate radiated EQ energy using NIED Hi-Net data
# Erg Calculations follow Boatwright et al. (2002)
# Japan regional correction values originate from Yoshimoto et al. (1993)
# Prints output : Station[1]  Channel[2]  Energy[3]   station_coord[4]   origin_coordinates[5]
# L.Barama November 2019
#------------------------------------------------------------------------------------------------
# Import Neccessary Modules
import sys
import glob
import os
from imp_bojoy import g2geomsprd
from imp_bojoy import atten
from imp_bojoy import imp_bojoy
from obspy import read
from obspy import read_inventory
from obspy.geodetics import locations2degrees
from obspy.taup import TauPyModel
import numpy as np
import math as m
from math import pi
import matplotlib.pyplot as plt
import scipy.integrate as integrate
from scipy.integrate import quad
from obspy.core import UTCDateTime
#------------------------------------------------------------------------
# How to Run:
#
#  the event number is comprised of EQ YEAR+MO_DAY+HR (UTC not JST!!)
#------------------------------------------------------------------------
##################################################################
# Set Constant values for Energy calculation:
# density and  S-wave velocity at receiver upper crustal values from (Ozel et al., 1999)
# Free-Surface Amplification (S) and root mean square S-wave radiation pattern coefficient (Fs)
# (from Boatwright et al. 2002, pg.1242)
S = 2
Fs = np.sqrt(2/5)
density = 2700    #  KG/m^3
beta = 3400       # S-wave velocity m/s
model = TauPyModel(model="iasp91")
#Define output path to output file with energy results
#outpath = "/home/lbarama/NIED/results/"
#outfile_name = Event+'.out'
#fout =  open(os.path.join(output_path,outfile_name), 'w') #Makes output files for each station
#list = open("Energy.out",'w')
# read traces in mseed files ( note these are velocity seismograms in nm/s)
file = "Events.edited"
with open(file, 'r') as read_file:
    for line in read_file:
        info = line.split()
        origlat = float(info[7])
        origlong = float(info[8])
        origdep = float(info[9]) # in kilometers
        origin_time = UTCDateTime(info[0]+'-'+info[1]+'-'+info[2]+' '+info[3]+':'+info[4]+':'+info[5])
        EventID = str(info[0]+info[1]+info[2]+info[3]+info[4])
        # Define the path to sac files and poles and zeros files /home/lbarama/NIED/Japan_data/1904110818
        path_to_sac = "/home/lbarama/NIED/Japan_data/"+EventID+"/sacfiles/"
        path_to_resp = "/home/lbarama/NIED/Japan_data/"+EventID+"/PolenZeros/"
        list = open("Energy."+EventID+".out",'w')
        for waveform in glob.glob(path_to_sac+'*SAC'):
            #print(waveform)
# Check that the waveform exists
            try:
                stream = read(waveform)
            except:
                continue
            net = stream[0].stats.network
            station = stream[0].stats.station
            loc = stream[0].stats.location
            channel = stream[0].stats.channel
            tr = stream[0]
            df = tr.stats.sampling_rate
            stalat = tr.stats.sac.stla #coordinates.latitude
            stalong = tr.stats.sac.stlo #coordinates.longitude
            start_time = tr.stats.starttime #start of trace
            # Find distance between station and event
            distance_degrees_orig_to_station  = locations2degrees(origlat,origlong,stalat,stalong)
            #Convert epicentral distance degrees to meters
            dist_m = 6371000.*(pi/180)*distance_degrees_orig_to_station
            # Find the predicted P-wave arrival time wrt the origin time (in seconds)
            predict_arrival = model.get_travel_times(source_depth_in_km = origdep, distance_in_degree = distance_degrees_orig_to_station, phase_list=["P"])
          # Need to get seconds (original output of the function above is the arrival time and the corresponding phase
            pred_Parrival_time = predict_arrival[0].time
            #print(pred_Parrival_time)
          # Re-define start of trace predicted P-arrival and 60 seconds after
            # Copy the traces! (To keep original data unchanged)
            tr_copy = tr.copy()
            tr_copy.trim(starttime=start_time+pred_Parrival_time, endtime=start_time+pred_Parrival_time+60)
            N = len(tr_copy)
            # Remove Mean and Take Bandpass filter 1 to 10 Hz
            tr_demean = tr_copy.detrend(type='demean')
            # could detrend as well
            v_filt = tr_demean.filter(type='bandpass',freqmin=1,freqmax=10, corners=4, zerophase = True) # zerophase = True, apply filter once forwards and once backwards.This results in twice the filter order but zero phase shift in the resulting filtered trace
            # Use HANN window Taper
            # To get the same results as the default taper in SAC, use max_percentage=0.05 and leave type as hann.
            v_taper = v_filt.taper(0.05,type='hann', max_length = None, side = 'both') # we're assuming the amount is per-side.
            v_costap = v_filt.taper(0.1,type='cosine', max_length = None, side = 'both')
            #convert velocity from nm/s to m/s
            v = v_taper.data*(1E-9)
            # Take fft and convert to angular frequency (mult by 2pi)
            # The abolute value of fft output gives the amplitude of the velocity spectrum (and only positive(realhalf))
            n = np.array(v).size
            fft_amplitude = abs((np.fft.fft(v)))[:n // 2]
            #The fftfreq function generates a list of frequencies that correspond to the components of the Fourier transform.
            # It gives $N$ values in the interval (-0.5,0.5).
            # To convert to the actual frequency, you need to divide by the sampling interval in time.
            n = np.array(v).size
            timestep = 1/df  #inverse of sampling rate
            # We want to include only the positive half of the frequencies to define frequencies freq
            freq = np.fft.fftfreq(n, d=timestep)
            freq_pos  = np.array(freq[:n // 2])
            n2 = np.array(v**2).size
            # We want to include only the positive half of the frequencies to define frequencies freq
            freq2 = np.fft.fftfreq(n, d=timestep)
            freq_pos2  = np.array(freq[:n // 2])
            #print(len(freq_pos))
            #print(min(freq),max(freq))
            #print(min(freq_pos),max(freq_pos))
            # Now Correct for Intrument Response!
            # Open and read from Corresponding PZ file
            #print(station+"."+channel+"."+"SAC_PZ")
            pzfile = open(path_to_resp+station+"."+channel+"."+"SAC_PZ","r")
            lines = pzfile.readlines()
            zerosline = lines[0]
            info = zerosline.split()
            nzero = int(info[1])
            polesline = lines[1]
            info = polesline.split()
            npole = int(info[1])
            polesline1 = lines[2]
            info = polesline1.split()
            preal1 = float(info[0])
            pcomp1 = float(info[1])
            polesline2 = lines[3]
            info = polesline2.split()
            preal2 = float(info[0])
            pcomp2 = float(info[1])
            #print(int(nzero),int(npole))
            #print(preal1,pcomp1)
            #print(preal2,pcomp2)
            #for waveform in glob.glob(polesnzeros+'/_PZ')
            Gain = 9.999484E-01
        #    for i in range(nzero):
        #        value = [0 + 0j,0 + 0j,0 + 0j,0 + 0j,0 + 0j]
        #        zeros = value[i]
            #print(zeros)
            zpole = complex(preal1, pcomp1),complex(preal2, pcomp2)
            f = freq_pos
            w = 2*pi*f
            '''
            for omega in range(len(w)):
                zw = complex(0,omega)
                for omega in range(len(w)):
                    zero = complex(0,omega)
                for i in range(0,nzero):
                    znumerator = 0
                    znumerator = znumerator*(zw -zero)
            #    for i in range(0,npole):
            #        pdenom = 0
            #        pdenom = pdenom*(zw-zpole[i])
            #zres = np.array(Gain*znumerator/pdenom)
            '''
            def paz_2_amplitude_value_of_freq_resp(paz, freq):
                """
                Returns Amplitude at one frequency for the given poles and zeros

                :param paz: Given poles and zeros
                :param freq: Given frequency

                The amplitude of the freq is estimated according to "Of Poles and
                Zeros", Frank Scherbaum, p 43.

                .. rubric:: Example

                >>> paz = {'poles': [-4.44 + 4.44j, -4.44 - 4.44j],
                ...        'zeros': [0 + 0j, 0 + 0j],
                ...        'gain': 0.4}
                >>> amp = paz_2_amplitude_value_of_freq_resp(paz, 1)
                >>> print(round(amp, 7))
                0.2830262
                """
                jw = complex(0, 2 * pi * freq)  # angular frequency
                fac = complex(1, 0)
                for zero in paz['zeros']:  # numerator
                    fac *= jw - zero
                for pole in paz['poles']:  # denominator
                    fac /= jw - pole
                return abs(fac) #* paz['gain']
                #Don't include Gain for Hinet data b/c Gain/Sensitivity is already removed when files are converted to sac

            paz = {'poles': zpole,'zeros': [0 + 0j, 0 + 0j],'gain': Gain}
            zres = []
            for f in range(len(freq_pos)):
                zres.append(paz_2_amplitude_value_of_freq_resp(paz, f))
            #print(len(np.array(zres)),np.array(zres))
            # Instrument Response Corrected velocity spectrum!!
            Vel_fft_amp = fft_amplitude/np.array(zres)
            #print(Vel_fft_amp)

            def id_nearest(array,scalar):      # reports index nearest value
                idx = (np.abs(array - scalar)).argmin()
                return idx

            f=freq_pos
            #print(np.array_equal(f,abs(freq_pos)))
            fbounds=np.array([2,10])
            fbids=np.zeros(fbounds.size)
            fbids[0]=id_nearest(f,fbounds[0])
            fbids[1]=id_nearest(f,fbounds[1])
            #print(fbids)
            f_cut=np.array(f[int(fbids[0]):int(fbids[1])])
            fft_amp_cut = np.array(Vel_fft_amp[int(fbids[0]):int(fbids[1])])
            # r = source-reciever distance, ro = radius of focal sphere from equation 5a Boatwright et al. 2002
            r = np.sqrt((dist_m)**2 + (origdep*1000)**2) # *1000 to convert to meters
            ro = 30000   # in meters and an estimate, check later
            # Calculate Corrected Velcocity spectrum! Equation 11 from Boatwright et al. 2002
            f=freq_pos
            fbounds=np.array([2,10])
            fbids=np.zeros(fbounds.size)
            fbids[0]=id_nearest(f,fbounds[0])
            fbids[1]=id_nearest(f,fbounds[1])
            f_cut=np.array(f[int(fbids[0]):int(fbids[1])])
            fft_amp_cut = np.array(Vel_fft_amp[int(fbids[0]):int(fbids[1])])
            ts = []
            v_corr = []
            for i in range(len(f_cut)):
                  k = 0.04 # seconds,generic rock site near surface impedence
                  w = 2*pi*f_cut[i]
                  Q = atten(f_cut[i]) # #atten_corr = atten(freq)check because it was written for frequency input but in equation is w
                  tstar = r/(beta*Q)
                  geo_corr = g2geomsprd(r,ro,f_cut[i])  #geo_corr = g2geomsprd(r,ro,freq)  ;  #imp_corr = imp_bojoy(w)
                  v_corr.append(np.sqrt((imp_bojoy(w)/(density*beta)))*np.exp((w*k*0.5) + ((w*tstar*.5)))*(geo_corr/(S*Fs))*fft_amp_cut[i])
            v2 = (np.array(v_corr)**2)# square the corrected velocity power spectrum
# Calculate the ENERGYYYYY
            df = f_cut[2]-f_cut[1]
            dw = 2*pi*df
            a = 0
            b = len(v2)
            integral = np.trapz(v2, dx=dw)  # computes integral using Trapezoidal Rule
            E = 4*density*beta*(Fs**2)*integral
# Write results to outfile
            list.write( str("%10.2e"%(E)) + '   '+ str(int(dist_m)/1000) + '   '+ str(origdep) + '   '+ str(stalat) + '   '+ str(stalong)+'   '+ str(origlat)+'   '+str(origlong) +'   ' +str(EventID) +'  '+str(station+'.'+channel) +'  '+ '\n' )
            print(str("%10.2e"%(E)) + '   '+ str(int(dist_m)/1000) + '   '+ str(origdep) + '   '+ str(stalat) + '   '+ str(stalong)+'   '+ str(origlat)+'   '+str(origlong) +'   ' +str(EventID) +'  '+str(station+'.'+channel) +'  '+ '\n')
        list.close()
