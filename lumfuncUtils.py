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
from dnds_lumfunc import *

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
                family=None,dump=None,verbose=False):

    """
    Do this for all bins simultaneously
    Mimic pn_integral.pro:
    ;+                                                                              
    ;                                                                               
    ; Ketron 11/2012                                                                
    ;                                                                               
    ; Integration of eq. (7).                                                       
    ;                                                                               
    ;                                                                               
    ;    /smax               /max_bin                                               
    ;    |     dS * dN/dS *  |        dS' exp(-(S-S')/2/rms^2)                      
    ;    /smin               /min_bin                                               
    ;                                                                               
    ; Analytic solution of second integral is error functions.                      
    ;                                                                               
    ;-  
    """

    C=alpha=Smin=Smax=beta=S0=gamma=S1=delta=S2=-99.0

    nlaws=int(0.5*len(paramsList)-1)

    C=params[paramsList.index('C')]
    Smin=params[paramsList.index('S0')]
    alpha=params[paramsList.index('a0')]
    if nlaws > 1:
        beta=params[paramsList.index('a1')]
        S0=params[paramsList.index('S1')]
    if nlaws > 2:
        gamma=params[paramsList.index('a2')]
        S1=params[paramsList.index('S2')]
    if nlaws > 3:
        delta=params[paramsList.index('a3')]
        S2=params[paramsList.index('S3')]

    iSmax=int([i for i in paramsList if i.startswith('S')][-1][-1])
    Smax=params[paramsList.index('S%i'%iSmax)]

    noise=params[paramsList.index('noise')]

    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    II2 = numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        #sqDeg2srr=1.0
        if nlaws == 1:
            II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C*1.0e-6*10**(-6.0*alpha),alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,sqDeg2srr*area),Smin,Smax)[0]
            II2[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2srr*area),Smin/1.0e6,Smax/1.0e6)[0]
        elif nlaws in [2,3,4]:
            II[ibin]=integrate.quad(lambda S:powerLawFuncErfsS(S,nlaws,C,alpha,-99.0,beta,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,S0/1.0e6,gamma,S1/1.0e6,delta,S2/1.0e6,noise/1.0e6,sqDeg2srr*area),Smin/1.0e6,Smax/1.0e6)[0]

    return II

#-------------------------------------------------------------------------------
