"""
Support functions for bayestack, bayestackClasses and lumfunc

Jonathan Zwart
August 2015

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

@profile
def riceIntegral(p,p0low,p0high,sigma):
    zlo = rice(p/sigma,scale=sigma).cdf(p0low)
    zhi = rice(p/sigma,scale=sigma).cdf(p0high)
    integral=(zlo-zhi)/2.0
    return integral

@profile
def F(p,p0,sigma):
    """
    This is just the rice distribution
    To obtain Rayleigh, set p0=0.0
    """
    return rice(p0/sigma,scale=sigma).pdf(p)

@profile
def GofP(p,sigma,p0low,p0high):
    """
    This is the integral of the Rice distribution,
    which is then a function of p
    P0 is a nuisance parameter (integrated over [p0low,p0high])
    """
    return integrate.quad(lambda P0: F(p,P0,sigma),p0lo,p0hi)[0]

#-------------------------------------------------------------------------------

gofp=numpy.vectorize(GofP)
f=numpy.vectorize(F)

#-------------------------------------------------------------------------------

# Polarization utility functions

def dNdP0(P,params=None,paramsList=None,inta=None):
    """
    Source count model
    Eventually this can incorporate powerLawFuncS, P0Dist etc.
    """

    C=alpha=Pmin=Pmax=beta=P0=gamma=P1=delta=P2=-99.0

    nlaws=int(0.5*len(paramsList)-1)

    C=params[paramsList.index('C')]
    Pmin=params[paramsList.index('S0')]
    alpha=params[paramsList.index('a0')]
    if nlaws > 1:
        beta=params[paramsList.index('a1')]
        P0=params[paramsList.index('S1')]
    if nlaws > 2:
        gamma=params[paramsList.index('a2')]
        P1=params[paramsList.index('S2')]
    if nlaws > 3:
        delta=params[paramsList.index('a3')]
        P2=params[paramsList.index('S3')]

    iPmax=int([i for i in paramsList if i.startswith('S')][-1][-1])
    Pmax=params[paramsList.index('S%i'%iPmax)]

    noise=params[paramsList.index('noise')]

    if Pmin < P < Pmax:
        if inta is not None:
            return float(inta(P))
        else:
            return powerLawFuncWrap(nlaws,P,C,alpha,-99.0,beta,\
                                       Pmin,Pmax,P0,gamma,P1,delta,P2,1.0)
    else:
        return 0.0

#-------------------------------------------------------------------------------

def P0Dist(b0Array,n0Array):
    """
    Return an interpolation object for mapping the dN/P0 distribution
    to an arbitrary value
    Eventually do analytically (given Q0, U0 distributions) --- if possible??
    
    Call as e.g.
    
    n0,b0,pp=plt.hist(P0,bins=bins,color='g',label='P0')
    inta=P0Dist(b0,n0)
    height=inta(pvalue)
    """

    return interp1d(medianArray(b0Array),n0Array,bounds_error=False,fill_value=0.0)

#-------------------------------------------------------------------------------

def rices(p0,sigma_QU,Pbinlow,Pbinhigh):
    """
    This is the Rice equivalent of the "erfs" variable of old
    """
    return rice(p0/sigma_QU,scale=sigma_QU).cdf(Pbinhigh) \
      - rice(p0/sigma_QU,scale=sigma_QU).cdf(Pbinlow)

#-------------------------------------------------------------------------------

def IP(dNdP0,params,paramsList,Pbinlow,Pbinhigh,inta=None):
    """
    Single integral version of I

    Double, inseparable integral:
    P (measured, cf S_m) runs from Pbinlow -> Pbinhigh
    P0 (underlying, cf S) runs from Pmin -> Pmax
    
    The integrand is the product of dN/dP0 (count model) and F (the Rician)

    The Rice functions are calculated analytically
    Gives identical results in some simple test cases
    Faster than I
    Works!
    """

    Pmin=params[paramsList.index('S0')]
    iPmax=int([i for i in paramsList if i.startswith('S')][-1][-1])
    Pmax=params[paramsList.index('S%i'%iPmax)]
    sigma_QU=params[paramsList.index('noise')]
    
    return integrate.quad(lambda p0:dNdP0(p0,params=params,paramsList=paramsList,\
             inta=inta)*rices(p0,sigma_QU,Pbinlow,Pbinhigh),Pmin,Pmax)[0]

#-------------------------------------------------------------------------------


@profile
def calculateP3(params,paramsList,bins=None,area=None,\
                family=None,dump=None,verbose=False,inta=None):

    """
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
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        #sqDeg2srr=1.0
        #WORKING HERE
        II[ibin]=sqDeg2srr*area*IP(dNdP0(params=params,paramsList=paramsList),\
            noise,Smin,Smax,bins[ibin],bins[ibin+1],inta=inta)
        #II[ibin]=integrate.quad(lambda S:powerLawFuncErfsS(S,nlaws,C,alpha,-99.0,beta,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,S0/1.0e6,gamma,S1/1.0e6,delta,S2/1.0e6,noise/1.0e6,sqDeg2srr*area),Smin/1.0e6,Smax/1.0e6)[0]

    return II

#-------------------------------------------------------------------------------

