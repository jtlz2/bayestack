"""
HI functions and evolution thereof

Jonathan Zwart
February 2016

"""

import os,sys
import importlib
import glob
import numpy
from math import pi,e,exp,log,log10,isinf,isnan
from scipy import integrate,stats
from scipy.interpolate import interp1d
from scipy.special import erf
import matplotlib.pyplot as plt
from profile_support import profile
from utils import sqDeg2sr,sqrtTwo,find_nearest,medianArray,\
                           interpol,buildCDF,Jy2muJy,interpola
from cosmocalc import cosmocalc

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

@profile
def erfss(S,Sbinlow,Sbinhigh,ssigma):
    """
    For a Gaussian noise model
    """
    return 0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

#-------------------------------------------------------------------------------

def get_Lbins(sbins,z,dl):
    """
    Units:
        sbins - muJy
        dl - Mpc
        Lbins - W Hz^-1
    """
    Lbins = [pi*4 *(s*1e-32)* (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1) for s in sbins]
    return Lbins #W/Hz

#-------------------------------------------------------------------------------

    
def get_sbins(fbins,z,dl):
    """
    Units:
        fbins - W Hz^-1 [luminosity bins]
        dl - Mpc
        sbins - ***Jy
    """
    sbins = fbins/(pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    return sbins*1e26 #Jy

#-------------------------------------------------------------------------------

@profile
def schechter(L,Lstar,alpha,norm):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha - no units
        * norm - phi_star Mpc^-3 mag^-1
    Outputs:
        Lbins -
        log10(phi) -
    """

    Lr = L/Lstar
    #print L, Lstar, alpha, norm
    phi = norm *(Lr**alpha) *numpy.exp(-Lr)/Lstar
    return [L],[numpy.log10(phi)]

#-------------------------------------------------------------------------------

@profile
def doublepowerlaw(L,Lstar,alpha1,alpha2,C):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha1 - no units
        * alpha2 - no units
        * C - phi_star Mpc^-3 mag^-1
    Outputs:
        Lbins -
        log10(phi) -
    """

    Lr = L/Lstar
    #print Lbins, Lstar, alpha1,alpha2, C
    phi = C/ ((Lr)**alpha1 + (Lr)**alpha2)
    return [L],[numpy.log10(phi)]

#-------------------------------------------------------------------------------

@profile
def dNdS_LF(S,z,dl,params=None,paramsList=None,inta=None,area=None,family=None):
    """
    Source count model
    S is in Jy
    Cosmology is set globally
    """

    Lmin=Lmax=Lnorm=Lstar=Lslope=Lslope2=Lzevol=-99.0

    if family in ['LFsch','LFdpl']:
        Lmin=params[paramsList.index('LMIN')]
        Lmax=params[paramsList.index('LMAX')]
        Lnorm=params[paramsList.index('LNORM')]
        Lstar=params[paramsList.index('LSTAR')]
        Lslope=params[paramsList.index('LSLOPE')]
        Lzevol=params[paramsList.index('LZEVOL')]
    if family=='LFdpl':
        Lslope2=params[paramsList.index('LSLOPE2')]

    #dl=cosmocalc(z,H0=Ho,WM=wm)['DL_Mpc']
    [Smin,Smax]=get_sbins([Lmin,Lmax],z,dl) # Jy

    if Smin < S < Smax:
        if inta is not None:
            return float(inta(S))
        else:
            L=get_Lbins([1.0e6*S],z,dl)
            if family=='LFsch':
                Lbins,log10phi=schechter(L[0],Lstar,Lslope,Lnorm)
            elif family=='LFdpl':
                Lbins,log10phi=doublepowerlaw(L[0],Lstar,Lslope,Lslope2,Lnorm)
            #print ',',L[0],10**log10phi[0],Lstar,Lslope,Lnorm
            phi=(10**log10phi[0])*(1.0+z)**Lzevol
            #if phi==0.0: phi=1.0e-99 # Hmm
            return phi
    else:
        return 0.0


#-------------------------------------------------------------------------------

@profile
def IL(dNdS_LF,redshift,dl,params,paramsList,Sbinlow,Sbinhigh,\
       inta=None,area=None,family=None):
    """
    The integrand is the product of dN/dS_LF (count model) and G (the Gaussian)
    Schechter: ['LNORM','LSTAR','LSLOPE','LMIN','LMAX','LZEVOL'] + ['noise']
    Double PL: ['LNORM','LSTAR','LSLOPE','LSLOPE2','LMIN','LMAX','LZEVOL'] + ['noise']
    """

    Lmin=params[paramsList.index('LMIN')]
    Lmax=params[paramsList.index('LMAX')]
    #dl=cosmocalc(redshift,H0=Ho,WM=wm)['DL_Mpc']
    [Smin,Smax]=1.0e6*get_sbins([Lmin,Lmax],redshift,dl)
    sigma=params[paramsList.index('noise')]
#    print sigma,Sbinlow,Sbinhigh,Smin,Smax
    return integrate.quad(lambda S:dNdS_LF(S,redshift,dl,params=params,paramsList=paramsList,\
             inta=inta,area=area,family=family)*erfss(S,sigma/1.0e6,Sbinlow/1.0e6,Sbinhigh/1.0e6),\
                                        Smin/1.0e6,Smax/1.0e6)[0]

#-------------------------------------------------------------------------------

@profile
def calculateL3(params,paramsList,bins=None,area=None,\
                family=None,dump=None,verbose=False,inta=None,redshift=None,dl=None):

    """
    For LF,
    function to calculate mock data for a given power law or interpolation object
    """

    #dl=cosmocalc(redshift,H0=Ho,WM=wm)['DL_Mpc']
    
    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        #sqDeg2srr=1.0
        #print area,sqDeg2sr
        II[ibin]=IL(dNdS_LF,redshift,dl,params,paramsList,bins[ibin],bins[ibin+1],\
                    inta=inta,area=sqDeg2srr*area,family=family)

    return II

#-------------------------------------------------------------------------------
