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

def polynomialFuncErfsS(S,modelFunc,Smin,Smax,Sbinlow,Sbinhigh,ssigma,area):
    """
    """

    if S < Smin or S > Smax:
        return 0.0

    erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))
    
    return erfs*modelFunc*area

#-------------------------------------------------------------------------------

@profile
def powerLawFuncS(S,C,alpha,Smin,Smax,area):
    """
    Ketron's equation (9)
    """

    if S < Smin or S > Smax:
        return 0.0
    else:
        n = C * (S ** alpha) * area
        return n

#-------------------------------------------------------------------------------


@profile
def powerLawFuncDoubleS(S,C,alpha,D,beta,Smin,Smax,S0,area):
    """
    """

    # Try a broken power law a la wikipedia:
    #http://en.wikipedia.org/wiki/Power_law#Broken_power_law

    if S <= Smin or S > Smax:
        return 0.0
    elif Smin < S <= S0:
        n = C * (S ** alpha) * area
    elif S0 < S <= Smax:
        n = C * S0**(alpha-beta) * (S ** beta) * area
    return n

#-------------------------------------------------------------------------------

@profile
def powerLawFuncTripleS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,area):

    if S <= Smin or S > Smax:
        return 0.0
    elif Smin < S <= S0:
        n = C * (S ** alpha) * area
    elif S0 < S <= S1:
        n = C * S0**(alpha-beta) * (S ** beta) * area
    elif S1 < S <= Smax:
        n = C  * S0**(alpha-beta) * S1**(beta-gamma) * (S ** gamma) * area

    return n

#-------------------------------------------------------------------------------

@profile
def powerLawFuncQuadS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,delta,S2,area):

    if S <= Smin or S > Smax:
        return 0.0
    elif Smin < S <= S0:
        n = C * (S ** alpha) * area
    elif S0 < S <= S1:
        n = C * S0**(alpha-beta) * (S ** beta) * area
    elif S1 < S <= S2:
        n = C  * S0**(alpha-beta) * S1**(beta-gamma) * (S ** gamma) * area
    elif S2 < S <= Smax:
        n = C  * S0**(alpha-beta) * S1**(beta-gamma) * \
                 S2**(gamma-delta)* (S**delta) * area

    return n

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
def calculateI(params,paramsList,bins=None,area=None,
                family=None,dump=None,verbose=False,model=None):

    """
    pn_integral, but for various different function families
    """

    if family=='ppl':
        return calculateI3(params,paramsList,bins=bins,area=area,\
                family=family,dump=dump,verbose=verbose)

    elif family=='poly':
        Smin=params[paramsList.index('S0')]
        Smax=params[paramsList.index('S1')]
        coeffs=[params[paramsList.index(p)] for p in paramsList if p.startswith('p')]
        ncoeffs=len(coeffs)
        noise=params[paramsList.index('noise')]

    for i in range(6): print params[i]
    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        print ibin,bins[ibin],bins[ibin+1],II[ibin],model.eval(bins[ibin],coeffs)
        II[ibin]=integrate.quad(lambda S:polynomialFuncErfsS(S,model.eval(S,coeffs),Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2srr*area),Smin/1.0e6,Smax/1.0e6)[0]

    return II

#-------------------------------------------------------------------------------

@profile
def calculateI3(params,paramsList,bins=None,area=None,\
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

    iSmax=int([i for i in paramsList if i[0]=='S'][-1][-1])
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
