"""
Support functions for bayestack, bayestackClasses and lumfunc
"""

import numpy
from math import pi,e,exp,log,log10,isinf,isnan
from scipy import integrate,stats
from scipy.interpolate import interp1d
from scipy.special import erf
from profile_support import profile
from utils import sqrtTwo,find_nearest,medianArray,interpol,buildCDF
from bayestack_settings import * # <-- generalize, localize

#-------------------------------------------------------------------------------

@profile
#def simulate(family,params,bins,nsources=None,noise=None,\
#             output='temp.txt',seed=None,dump=False,\
#             verbose=False,output=None,version=2):
def simulate(seed=None,N=None,noise=None,output=None,dump=None,\
             version=2,verbose=False):
    """
    Based on lumfunc.simtable()
    Specify family + parameters
    Specify number of sources
    Build CDF (or set up function)
    Draw deviates
    Sample CDF given deviates
    Add noise (None or some value)
    Bin
    Write out
    Return

    Look at simulate.ipynb for an example run
    """

    if version < 2: return '***Unsupported!!'
    
    if seed is not None:
        numpy.random.seed(seed=SEED_SIM)

    function = lambda S:S**2

    # Set up the 'rough' array
    Smin=0.0 # uJy
    Smax=100.0 # uJy
    gridlength=10000
    Ss=numpy.linspace(Smin,Smax,gridlength)
    values=numpy.array([function(ix) for ix in Ss])

    # Build the CDF
    CDF=buildCDF(values)
    # Create the interpolant object
    sampler=interp1d(CDF,Ss)

    # Draw the random deviates
    R = numpy.random.rand(N)
    F=sampler(R)

    # Integrate the original function
    A=integrate.quad(function,Smin,Smax)[0]

    # Bin the random samples
    bins=numpy.linspace(Smin,Smax,40)
    dbin=bins[-1]-bins[-2]
    E=numpy.histogram(F,bins=bins)[0]
    # And calculate their area
    C=integrate.trapz(E)*dbin

    # Gunpowder, treason and....
    #plt.xlim(0.0,100.0)
    #plt.xlabel('S / $\mu$Jy')
    #plt.hist(F,bins=bins)
    #plt.plot(Ss_fine,values*C/A,'r')

    # Dump noiseless fluxes to file
    if dump is not None:
        dumpf=dump
        numpy.savetxt(dumpf,F)
        print 'Draws (noiseless) are in %s' % dumpf

    # Add noise if requested
    if noise is not None:
        numpy.random.seed(seed=SEED_SIM)
        F+=numpy.random.normal(0.0,noise,N)

    # Dump noisy fluxes to file
    if dump is not None:
        noisydumpf='%s_noisy.txt' % dumpf.split('.')[0]
        numpy.savetxt(noisydumpf,F)
        print 'Draws (noisy) are in %s' % noisydumpf
        print 'Minimum flux in catalogue = %f' % F.min()

    # Bin up the fluxes
    counts=numpy.histogram(R,bins=bins)[0]
    print '-> %i/%i objects observed in total (bin scattering)\n' % (counts.sum(),N)

    # Calculate differential counts
    idl_style=False
    dn_by_ds,dn_by_ds_eucl,dn_by_ds_errs,dn_by_ds_b,dn_by_ds_b_errs=\
      calculateDnByDs(1.0e-6*bins,counts,idl_style=idl_style,return_all=True)

    median_bins=medianArray(bins) # uJy
    if output is not None:
        outputf=output
        #data=numpy.zeros((len(bins),3))
        s=open(outputf,'w')
        if version < 2:
            header='# bin_median ksRaw ksNoisy'
        else:
            header='# bin_low_uJy bin_high_uJy bin_median_uJy Ns_tot_obs dnds_srm1Jym1 dnds_eucl_srm1Jy1p5 delta_dnds_eucl_lower_srm1Jy1p5 delta_dnds_eucl_upper_srm1Jy1p5 corr Ncgts_degm2 dNcgts_lower_degm2 dNcgts_upper_degm2'
        s.write('%s\n'%header)
        if verbose: print header
        for ibin in range(nbins-1):
            if version < 2:
                line='%f %i %i' % (median_bins[ibin],-99.0,counts[ibin])
            else:
                # See binner.py for dependence on BIN_CORRS/CORR_BINS (omitted here)
                line='%f %f %f %i %e %f %f %f %f %i %i %i' % \
                  (bins[ibin],bins[ibin+1],median_bins[ibin],round(counts[ibin]),\
                   dn_by_ds[ibin]/(sqDeg2sr*AREA_SIM),\
                   dn_by_ds_eucl[ibin]/(sqDeg2sr*AREA_SIM),\
                   dn_by_ds_errs[ibin]/(sqDeg2sr*AREA_SIM),\
                   dn_by_ds_errs[ibin]/(sqDeg2sr*AREA_SIM),\
                   1.00,\
                   round(counts[ibin:].sum()*1.00/AREA_SIM),\
                   round(numpy.sqrt(counts[ibin:].sum()*1.00/AREA_SIM)),\
                   round(numpy.sqrt(counts[ibin:].sum()*1.00/AREA_SIM)))
            s.write('%s\n'%line)
            if verbose: print line
        print counts.sum()
        s.close()
    
    return F

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
def polynomialFuncErfsS(S,S_1,coeffs,Smin,Smax,Sbinlow,Sbinhigh,ssigma,aarea):
    """
    """

    if S < Smin or S > Smax:
        return 0.0

    erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

    return erfs * polyFunc(S,S_1,coeffs) * aarea

#-------------------------------------------------------------------------------

@profile
def polesFuncErfsS(S,pole_posns,coeffs,Sbinlow,Sbinhigh,ssigma,aarea):
    """
    """

    if S < pole_posns[0] or S > pole_posns[-1]:
        return 0.0

    erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

    return erfs * polesFunc(S,pole_posns,coeffs) * aarea

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
def polyFunc(S,S_1,c):
    """
    """
    exponent=0.0
    for n in range(len(c)):
        exponent += c[n] * (numpy.log10(S/S_1)**n)

    return 10**exponent

#-------------------------------------------------------------------------------

@profile
def polesFunc(S,pole_posns,coeffs):
    """
    """

    
    #print coeffs
    #print S,idx,Snearest#,pole_posns[0],pole_posns[-1]
    #sys.exit(0)
    #print idx,pole_posns[idx],coeffs[idx]
#    counts=calculateDnByDs(pole_posns,coeffs[:-1],eucl=False,verbose=False,\
#                    idl_style=True,errors=False,bright=False,bright_errors=False,
#                    return_all=False)

    idx,Snearest=find_nearest(pole_posns[:-1],S)
    #if idx >= len(coeffs)-1: return counts[idx-1]
    #return counts[idx]
    return coeffs[idx]/(pole_posns[idx+1]-pole_posns[idx])

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
#        func=numpy.poly1d(list(reversed(coeffs)))
#        ncoeffs=len(coeffs)
        S_1=1.0 # ref flux
        noise=params[paramsList.index('noise')]

        #for i in range(6): print params[i]
        nbins=len(bins)-1
        II = numpy.zeros(nbins)
        for ibin in xrange(nbins):
            II[ibin]=integrate.quad(lambda S:polynomialFuncErfsS(S,S_1,coeffs,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2sr*area),Smin/1.0e6,Smax/1.0e6)[0]
            #print ibin,bins[ibin],bins[ibin+1],II[ibin]
        return II

    elif family=='bins':
        coeffs=[params[paramsList.index(p)] for p in paramsList if p.startswith('b')]
        noise=params[paramsList.index('noise')]
        pole_posns=numpy.logspace(-1,numpy.log10(85.0),len(coeffs)+1)
        #pole_posns=numpy.linspace(1.0,85.0,11)
        assert(len(coeffs)==len(pole_posns)-1), '***Mismatch in number of poles!!'

        nbins=len(bins)-1
        II = numpy.zeros(nbins)
        for ibin in xrange(nbins):
            II[ibin]=integrate.quad(lambda S:polesFuncErfsS(S,pole_posns/1.0e6,coeffs,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2sr*area),pole_posns[0]/1.0e6,pole_posns[-1]/1.0e6)[0]
            #print II[ibin]
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
