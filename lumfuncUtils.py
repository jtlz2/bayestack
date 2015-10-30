"""
Support functions for bayestack, bayestackClasses

Luminosity functions and evolution thereof

Jonathan Zwart
October 2015

"""

import os,sys
import importlib
import glob
import numpy
from math import pi,e,exp,log,log10,isinf,isnan
from scipy import integrate,stats
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.stats import rice,rayleigh
import matplotlib.pyplot as plt
from profile_support import profile
from utils import sqDeg2sr,sqrtTwo,find_nearest,medianArray,\
                           interpol,buildCDF,Jy2muJy,interpola
from countUtils import powerLawFuncWrap
from cosmocalc import cosmocalc
import dnds_lumfunc

if 'chains' in sys.argv[-1]:
    potential_settings=glob.glob(os.path.join(sys.argv[-1],'*settings*py'))
    assert(len(potential_settings)==1), '***More than one potential settings file!'
    settingsf='.'.join([sys.argv[-1],potential_settings[0].split('/')[-1].split('.')[-2]])
else:
    settingsf=sys.argv[-1].split('.')[-2]

print '%s is using %s' % (__name__,settingsf)
try:
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)
except:
    print '***Warning: Settings not loaded'

#-------------------------------------------------------------------------------

# LF utility functions
def erfss(S,Sbinlow,Sbinhigh,ssigma):
    return 0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

#-------------------------------------------------------------------------------

@profile
def dNdS_LF(S,z,params=None,paramsList=None,inta=None,area=None):
    """
    Source count model
    S is in Jy
    """

    Lmin=Lmax=Lnorm=Lstar=Lslope=Lzevol=-99.0
    Lmin=params[paramsList.index('LMIN')]
    Lmax=params[paramsList.index('LMAX')]
    Lnorm=params[paramsList.index('LNORM')]
    Lstar=params[paramsList.index('LSTAR')]
    Lslope=params[paramsList.index('LSLOPE')]
    Lzevol=params[paramsList.index('LZEVOL')]

    dl=cosmocalc(z,H0=Ho,WM=wm)['DL_Mpc']
    [Smin,Smax]=dnds_lumfunc.get_sbins([Lmin,Lmax],z,dl)

    if Smin < S < Smax:
        if inta is not None:
            return float(inta(S))
        else:
            L=dnds_lumfunc.get_Lbins([S],z,dl)
            Lbins,log10(phi)=schechter(L,Lstar,Lslope,Lnorm)
            return (1.0+z)**Lzevol*10**phi[0]
    else:
        return 0.0


#-------------------------------------------------------------------------------

@profile
def calculateIL(params,paramsList,bins=None,area=None,
                family=None,dump=None,verbose=False,model=None):

    """
    pn_integral, but for luminosity function families
    Flux arguments to this function are in uJy
    ['LNORM','LSTAR','LSLOPE','LMIN','LMAX','LZEVOL'] + ['noise']
    """

    if family=='LFsch':
        Lmin=params[paramsList.index('LMIN')]
        Lmax=params[paramsList.index('LMAX')]
        Lnorm=params[paramsList.index('LNORM')]
        Lstar=params[paramsList.index('LSTAR')]
        Lslope=params[paramsList.index('LSLOPE')]
        Lzevol=params[paramsList.index('LZEVOL')]
        noise=params[paramsList.index('noise')]
        redshift=?????
        #for i in range(6): print params[i]
        nbins=len(bins)-1
        II = numpy.zeros(nbins)
        for ibin in xrange(nbins):
            II[ibin]=integrate.quad(lambda S:polynomialFuncErfsS(S,S_1,coeffs,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2sr*area),Smin/1.0e6,Smax/1.0e6)[0]
            #print ibin,bins[ibin],bins[ibin+1],II[ibin]
        return II


#-------------------------------------------------------------------------------

@profile
def calculateL3(params,paramsList,bins=None,area=None,\
                family=None,dump=None,verbose=False,inta=None,doRayleigh=False):

    """
    For LF,
    function to calculate mock data for a given power law or interpolation object
    """

    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        #sqDeg2srr=1.0
        #print area,sqDeg2sr
        II[ibin]=IL(dNdS_LF,params,paramsList,bins[ibin],bins[ibin+1],inta=inta,area=sqDeg2srr*area)

    return II

#-------------------------------------------------------------------------------
