"""
Support functions for bayestack, bayestackClasses and lumfunc
"""

import numpy
from math import pi,e,exp,log,log10,isinf,isnan
from scipy import integrate,stats
from scipy.special import erf
from profile_support import profile
from utils import sqrtTwo
from bayestack_settings import * # <-- generalize, localize

#-------------------------------------------------------------------------------

@profile
def calculateDnByDs(bins,counts,eucl=False,verbose=False,idl_style=True,
                    errors=False,bright=False,bright_errors=False,
                    return_all=False):
    """

    This function expects bins to be in uJy
    The output units are uJy^-1, or uJy^1.5 if eucl=True

    The function is area-agnostic: sr^-1 -> sr^-1, deg^-2 -> deg^-2 etc.

    Be careful with units:
       e.g. if bins are in uJy, dn_by_ds will be uJy^-1   = 10^6 Jy
                                dn_ny_ds_eucl in uJy^1.5  = 10^-9 Jy
    etc.
    """

    dn=counts
    Smed=medianArray(bins)
    ds=numpy.absolute(numpy.gradient(Smed))

    if idl_style:
        # i.e. the way DERIV.PRO does it
        ds[0] = abs(-3.0*Smed[0] + 4.0*Smed[1] - Smed[2] / 2.0)
        ds[-1] = abs(3.0*Smed[-1] - 4.0*Smed[-2] + Smed[-3] / 2.0)

    if verbose:
        print 'dn',dn
        print 'ds',ds

    if return_all:
        return dn/ds, (Smed**2.5) * dn/ds,\
          (Smed**2.5) * numpy.sqrt(dn)/ds,\
          (Smed**2.0) * dn/ds,\
          (Smed**2.0) * numpy.sqrt(dn)/ds

    if eucl:
        dn_by_ds = (Smed**2.5) * dn/ds
    elif errors:
        dn_by_ds = (Smed**2.5) * numpy.sqrt(dn)/ds
    elif bright:
        dn_by_ds = (Smed**2.0) * dn/ds
    elif bright_errors:
        dn_by_ds = (Smed**2.0) * numpy.sqrt(dn)/ds
    else:
        dn_by_ds = dn/ds

    return dn_by_ds

#-------------------------------------------------------------------------------


@profile
def powerLawFuncErfsS(S,nlaws,C,alpha,D,beta,Smin,Smax,\
                      Sbinlow,Sbinhigh,S0,gamma,S1,delta,S2,ssigma,area):
    """
    Ketron's equation (5) + (9)

    ; Integration of eq. (7).                                                       
    ;
    ;    /smax               /max_bin                                               
    ;    |     dS * dN/dS *  |        dS' exp(-(S-S')**2/2/rms^2)                      
    ;    /smin               /min_bin                                               
    ;                                                                               
    ; Analytic solution of second integral is error functions.                      
    ;
    """

    if S < Smin or S > Smax:
        return 0.0

    erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))
    return erfs * powerLawFuncWrap(nlaws,S,C,alpha,D,beta,\
                                       Smin,Smax,S0,gamma,S1,delta,S2,area)

#-------------------------------------------------------------------------------

@profile
def powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Smin,Smax,S0,gamma,S1,delta,S2,area):

    if nlaws == 1:
        return powerLawFuncS(S,C,alpha,Smin,Smax,area)
    elif nlaws == 2:
        return powerLawFuncDoubleS(S,C,alpha,-99.0,beta,Smin,Smax,S0,area)
    elif nlaws == 3:
        return powerLawFuncTripleS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,area)
    elif nlaws == 4:
        return powerLawFuncQuadS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,delta,S2,area)

    return

#-------------------------------------------------------------------------------

@profile
def powerLawFuncErrorFn(Si,C,alpha,Smin,Smax,Sbinlow,Sbinhigh,noise,area):

    """
    Units:
            Si, Sbinlow, Sbinhigh, noise - must be the same as each
                                             other, i.e. Jy or uJy
            erfs, alpha - dimensionless
            C - /([Si]*[area]), e.g. Jy^-1 sr^-1, uJy^-1 deg^-2, or some mix
            Hence units of C are tied to units of Si
            area is often just set to 1.0 anyway
    """

    if Si < Smin or Si > Smax:
        return 0.0

    erfs = erf((Si-Sbinlow)/(sqrtTwo*noise)) - erf((Si-Sbinhigh)/(sqrtTwo*noise))
    erfs *= 0.5
    n = C * Si**alpha * erfs * area

    return n

#-------------------------------------------------------------------------------

@profile
def calculateI3(C,alpha,Smin,Smax,area,noise=None,dump=None,\
                verbose=False,nlaws=1,D=None,beta=None,S0=None,\
                gamma=-99.0,S1=-99.0,delta=-99.0,S2=-99.0,bins=None):
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

    nbins=len(bins)
    
    II = numpy.zeros(nbins-1)
    II2 = numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        #sqDeg2srr=1.0
        if nlaws == 1:
            II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C*1.0e-6*10**(-6.0*alpha),alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,sqDeg2srr*SURVEY_AREA),Smin,Smax)[0]
            II2[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2srr*SURVEY_AREA),Smin/1.0e6,Smax/1.0e6)[0]
        elif nlaws in [2,3,4]:
            II[ibin]=integrate.quad(lambda S:powerLawFuncErfsS(S,nlaws,C,alpha,D,beta,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,S0/1.0e6,gamma,S1/1.0e6,delta,S2/1.0e6,noise/1.0e6,sqDeg2srr*SURVEY_AREA),Smin/1.0e6,Smax/1.0e6)[0]

    return II

#-------------------------------------------------------------------------------
