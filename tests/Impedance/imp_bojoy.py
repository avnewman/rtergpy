#!/usr/bin/env python3

def g2geomsprd(r,r0,f): #this is the g(r/r0) spreading parameter
    """
    at regional distances, waveforms are highly dispersed, with
    body waves dominating at short periods/high-f (f>2Hz)
    Lg-waves dominate at intermediate (0.3≤f≤1Hz)
    Surface at long period/low-f(<0.2Hz)
    """
    if r <=r0:
        g2=r*r
    else:
        if f<=0.2:   # surface wave
            gamma=0.5
        elif f<=0.25:  # Lg waves
            gamma=0.5+2*log(5*f)
        else:
            gamma=0.7  # body waves
        g2=r0*(r/r0)**(2*gamma)
    return g2

def atten(f): # This is the intrinsic/anelastic attenuation only.  Scattering attenuation will be measured as part of the total energy
    #if f<=1.5:
    #    bigQ=400
    #else:
    #    bigQ=400*(f/1.5)**0.6
    # Yoshimoto et al 1993
    bigQ = 1/(0.012*f**-0.73)  # they report in Qs^-1  # valid only around 1 - 24 Hz
    return bigQ

def imp_bojoy(w):
    """
    Calculate Impedence following Boatwright and Choy 2002's adaptation (Equation 9)
    of Boore and Joyner (1997). This algorithm will calculate site-side Impedence using
    generic "hard rock" density and shear parameters integrated down to 8 km depth.
    Depth of penetration (z) is assumed to be will be determined by frequency of signal,
    such that the wave will only sample the top quarter of the period (T/4).

    Inputs:
      x = Distance (m)
      w = angular frequency (2 pi f)
    Output:
      imp = Impedence (kg/m^2/s)
    """

    from math import pi
    import numpy as np
    from numpy import array
    from scipy.interpolate import interp1d
    import sys

#  depth constraints for allowing rock parameterization
    z0=0
    zmax=30.0  # [km] (currently going way beyond our limits) the maximum depth generic rock descp = 8 km
    zmax=8.0   # [km] this is the max depth of the generic rock description
    nsteps=500 # number of points to predefine parameters with. should be 100-500. Omega will be picked from
               # an iteration.

    def genrockbeta(z):
        from numpy import array
        beta=np.zeros(z.size)
        Beta=np.array([
            [   0.,  .001,  .245,    0],
            [ .001,   .03, 2.206, .272],
            [ .030,  .190, 3.542, .407],
            [ .190, 4.000, 2.505, .199],
            [4.000, 8.000, 2.927, .086]])  #// Table 1 from Boore and Joyner 1997 (Generic Rock Shear velocity)
           # zmin, zmax, beta_0, ^exponent:   solution in km/s
        i=0
        for zi in np.nditer(z):
            if zi <= Beta[0,1]:
                beta[i]=Beta[0,2]*zi**Beta[0,3]
            elif zi <= Beta [1,1]:
                beta[i]=Beta[1,2]*zi**Beta[1,3]
            elif zi <= Beta [2,1]:
                beta[i]=Beta[2,2]*zi**Beta[2,3]
            elif zi <= Beta [3,1]:
                beta[i]=Beta[3,2]*zi**Beta[3,3]
            elif zi <= Beta [4,1]:
                beta[i]=Beta[4,2]*zi**Beta[4,3]
            else:
                beta[i]=Beta[4,2]*zi**Beta[4,3]  #extrapolate just don't use for very long periods
                #sys.exit("Frequency out of range.")
            i += 1  # pythons way
        return beta

    def id_nearest(array,value):      # reports index nearest value from array
        idx = (np.abs(array - value)).argmin()
        return idx




    # a number of numpy arrays to be populated with generic rock informaiton
    z=np.linspace(z0,zmax,nsteps,endpoint=True)  # create a priori linear distribution.
    beta=np.zeros(z.size)
    tp=np.zeros(z.size)
    omega=np.zeros(z.size)
    f=np.zeros(z.size)
    t_over4=np.zeros(z.size)
    density=np.zeros(z.size)

    beta=genrockbeta(z)  # call local program to calculate velocity profile for range of depths
    tp=z/beta           # period for those ranges
    #t_over4=tp/4.       # quarter wavelength window. (quarter period in Boatwright/Boore lingo)
    t_over4=pi/2/w      # Integrate from Z @ T/4.
    f[0]=1.e6           #  fix NaN
    f=1./tp             # Frequency (have to put that factor of 4 back in..)
    omega=2.*pi*f       # angular frequency
    density=(2.5+(beta-0.3)*((2.8-2.5)/(3.5-0.3)))*1000. # ugly equation from Boore and Joyner
    #print(z,beta,t_over4,tp,f,omega)  # (worth freq/z relationzhip)
    idx=id_nearest(omega,w)
    dz=np.diff(z) # is constant (or better be)
    #print(idx,z[idx],f[idx],omega[idx],t_over4[idx],density[idx])

    if isinstance(w, (int, float)):
        idx=id_nearest(omega,w*4)  # find position at T/4 (T=2pi/omega)
        imp = z[idx]*(np.sum(dz[:idx]/(density[:idx]*beta[:idx]))**-1) #
        #print(idx,z[idx],f[idx],omega[idx],t_over4[idx],density[idx], imp)

    elif len(w) > 1:
        imp=np.zeros(w.size)
        i=0
        for wi in np.nditer(w):
            idx=id_nearest(omega,wi*4) # find position at T/4 (T=2pi/omega)
            imp[i]=z[idx]*(np.sum(dz[:idx]/(density[:idx]*beta[:idx]))**-1)
            i+=1

    return imp
