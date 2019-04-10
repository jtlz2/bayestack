#!/usr/bin/env python
##!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
##!/usr/bin/python

"""
This is lumfunc.py
Jonathan Zwart
January 2014
"""

from mpi4py import MPI
import dill
MPI._p_pickle.dumps = dill.dumps
MPI._p_pickle.loads = dill.loads

import os,sys,shutil
import importlib
import ctypes
import time,threading
import numpy
import pymultinest
#if 'baltasar' not in os.getenv('PBS_O_HOST'):
#from matplotlib import pyplot as plt
from math import pi,e,exp,log,log10,isinf,isnan
from scipy import integrate,stats
from scipy.interpolate import interp1d
# Should be faster than math.erf
# -- see http://stackoverflow.com/questions/457408/is-there-an-easily-available-implementation-of-erf-for-python
from scipy.special import erf
from priors import Priors
from profile_support import profile
from utils import *
import pylab
import stacker



# See
# http://stackoverflow.com/questions/15768136/python-bulk-promote-variables-to-parent-scope
#param_file='settings.py'
# This needs to be here so that the functions can pick up the settings
# variables

if __name__ == '__main__' and len(sys.argv) < 2:
    print 'usage:'
    #print '       ./lumfunc.py NOISY  SIGMA  AREA  MODEL  LOUD  TRUTH  SEED NLAWS'
    #print '                    [True|False] [uJy] [sq. deg.]  [0|1|2] [True|False]'
    #print '                    [True|False] (suppressible by LOUD) [1|2]'
    #print '                    [-1 for clock]'
    print
    print 'with MPI [NB mpiEXEC on marchmain]:'
    print '       mpiexec -n NPROC ./lumfunc.py settings.py'
    print
    print 'Models:'
    print '       0  vary no parameters'
    print '       1  vary all 4 parameters'
    print '       2  vary C and alpha only'
    print '       3  vary Smin and Smax only'
    print
    sys.exit(0)

# This is to deal with the fact that __main__ is being reset by importlib(?)
# This section is EXTREMELY messy and to be sorted out eventually....
# Need to remove some circularities

__name_cached=__name__
param_file=sys.argv[-1]

#setf='%s' % param_file.split('.')[-2]
# Import the settings variables
if __name__=='__main__':
    setf='settings'
else:
    setf='%s'%param_file.split('.')[0]

# Plot context
if 'plot' in sys.argv[0]:
    setf='%s.settings'%param_file

# command line context
if __name__=='lumfunc':
    context='cmdline'
    setf='settings'

# If not running on cluster:
if MPI.COMM_WORLD.size != 0:
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)
__name__=__name_cached

#try:
#    execfile(param_file)
#except IOError:
#    from settings import *

#-------------------------------------------------------------------------------


#@profile
#def calculateSigmaCStar(synthBeamSolidAngleSr,fluxes,counts):
#    """
#    """
#    mask1=(fluxes<=Smax)
#    SigSqConf=synthBeamSolidAngleSr*integrate.trapz(fluxes[mask]**2*counts[mask],fluxes[mask])
#    return
#-------------------------------------------------------------------------------


@profile
def powerLaw2Array(nlaws,C,alpha,beta,S0,gamma,S1,delta,S2,Smin,Smax):

    """
    Function to aid testing of calculateConfusionNoiseSqArray()
    """

    fluxes=numpy.logspace(log10(Smin),log10(Smax),1000)
    f=lambda S:powerLawFuncWrap(nlaws,S,C,alpha,-99.0,beta,\
                Smin,Smax,S0,gamma,S1,delta,S2,SURVEY_AREA)
    counts=numpy.array(map(f,fluxes))

    return fluxes,counts

#-------------------------------------------------------------------------------

@profile
def calculateConfusionNoiseSqArray(synthBeamSolidAngleSr,fluxes,counts,Smin,Smax):

    """
    Calculate confusion noise given an array of S,dN/dS, cf Vernstrom
    """
    mask1=(fluxes<=Smax)
    #mask2=(fluxes>=Smin)
    mask=mask1
    #print fluxes[mask][0],fluxes[mask][-1]
    SigSqConf=synthBeamSolidAngleSr*integrate.trapz(fluxes[mask]**2*counts[mask],fluxes[mask])

    return SigSqConf

#-------------------------------------------------------------------------------

@profile
def calculateConfusionNoiseSq(synthBeamSolidAngleSr,nlaws,\
                            C,alpha,beta,S0,gamma,S1,Smin,Smax):

    """
    For a given source count model, calculate the squared confusion noise
    See JZ thesis eqn
    """

    DUMMY=1.0
    I=integrate.quad(lambda S:
                S**2*powerLawFuncWrap(nlaws,S,C,alpha,-99.0,beta,\
                                              Smin,Smax,S0,gamma,S1,DUMMY),\
                                              Smin,Smax)[0]

    SigSqConf = synthBeamSolidAngleSr * I

    return SigSqConf

#-------------------------------------------------------------------------------

@profile
def simBins(bins,nlaws,C,alpha,D,beta,Smin,Smax,S0,gamma,S1,area,\
            intAlgorithm='quad',verbose=False):
    """
    Simulate counts in bins given a defined power law
    Do this for ALL the bins
    THIS IS NOISELESS
    intAlgorithm='quad' or 'brute' (quad is faster and probably more accurate)
    NOT YET WORKING FOR AREAS != 10.0 SQ. DEG.
    Call as e.g.
        II=lumfunc.simBins(lumfunc.bins,1,40.0,-1.5,-99.0,-99.0,1.0,20.0,-99.0,10.0)
    """

    II = numpy.zeros(nbins-1)
    #bin_medians=medianArray(bins)

    resolution = 1000
    #print '# binlow binhigh binmedian ksRaw I diffI'
    for ibin in range(nbins-1):
        Sbinlow=bins[ibin]
        Sbinhigh=bins[ibin+1]
        I=0.0
        if intAlgorithm == 'brute':
            ds = (Sbinhigh - Sbinlow) / resolution
            for istep in xrange(resolution):
                Si = Sbinlow + (istep - 0.5) * ds
                # An extra factor of resolution seems to be needed here (cf Ketron)
                #I += powerLawFuncS(Si,C,alpha,Smin,Smax,area) * ds * resolution
                I += powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Smin,Smax,S0,gamma,S1,area)\
                    * ds * resolution

        elif intAlgorithm == 'quad':
            #I=integrate.quad(lambda S: powerLawFuncS(S,C,alpha,Smin,Smax,area),\
            #                 Sbinlow,Sbinhigh)[0]
            I=integrate.quad(lambda S:
                             powerLawFuncWrap(nlaws,S,C,alpha,D,beta,\
                                              Smin,Smax,S0,gamma,S1,area),\
                                              Sbinlow,Sbinhigh)[0]

            # WHY?????????
            I *= resolution

        #print Sbinlow,Sbinhigh,bin_medians[ibin],ksRaw[ibin],int(I),int(I-ksRaw[ibin])
        #print int(J),int(J-ksRaw[ibin])
        II[ibin] = I

    dn_by_ds=calculateDnByDs(bins,II)
    dn_by_ds_eucl=calculateDnByDs(bins,II,eucl=True)
    dn_by_ds_eucl_ksNoisy=calculateDnByDs(bins,ksNoisy,eucl=True)
    dn_by_ds_eucl_ksRaw=calculateDnByDs(bins,ksRaw,eucl=True)
    if verbose:
        print '# binlow binhigh binmedian ksRaw ksNoisy I diffI dnds dnds_eucl log10_dnds_eucl log10_dnds_eucl_ksNoisy log10_dnds_eucl_ksRaw'
        for ibin in range(nbins-1):
            print bins[ibin],bins[ibin+1],bin_medians[ibin],ksRaw[ibin],ksNoisy[ibin],
            print int(II[ibin]),int(II[ibin]-ksRaw[ibin]),
            print int(dn_by_ds[ibin]),int(dn_by_ds_eucl[ibin]),
            print log10(sqDeg2sr*dn_by_ds_eucl[ibin]),
            print log10(sqDeg2sr*dn_by_ds_eucl_ksNoisy[ibin])
            print log10(sqDeg2sr*dn_by_ds_eucl_ksRaw[ibin])

    return II

#-------------------------------------------------------------------------------

#bin_medians=medianArray(bins)

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
def C2N(C,alpha,Smin,Smax,area):

    """
    Integrate the power law and return the number of objects under it
    """

    N=int(round(integrate.quad(lambda S:powerLawFuncS(\
            S,C,alpha,Smin,Smax,area),Smin,Smax)[0]))

    return N
    
#-------------------------------------------------------------------------------

@profile
def sampleSKADS(skadsf='heywood/wilman_counts_1p4GHz_jz.txt'):
    """
    """

    skads=numpy.genfromtxt(skadsf)
    S=skads[:,0] # S
    dnds=numpy.power(skads[:,0],-2.5)*skads[:,3] # dnds
    N=skads[:,2] # N

    Smin=S[0]; Smax=S[-1]
    Smin=0.01e-6; Smax=85.0e-6
    iSmin,Smin=find_nearest(S,Smin)
    iSmax,Smax=find_nearest(S,Smax)
    #iSmin=0
    S=S[iSmin:iSmax]
    N=N[iSmin:iSmax]

    Np=numpy.zeros(len(N))
    #for i,n in enumerate(N):
    #    Np[i]=numpy.random.poisson(lam=N[i])
    #    print i,S[i],N[i],Np[i]

    print '%f -> %f' % (N.sum(),Np.sum())


    #np=numpy.zeros(len(N))
    #np = Np * (radioSynthOmegaSr) # Flux per pixel
    #print np
    radioSynthSqDeg=radioSynthOmegaSr/sqDeg2sr
    #N_j=(value from SKADS)*(area of pixel in square degrees)
    n=N*1.7361111111111115e-07#radioSynthSqDeg
    #n=N*radioSynthSqDeg/64.0
    #n=N*1.7361111111111115e-07/64.0
    npix=72000
    Sp=numpy.zeros(npix)
    nprunning=numpy.zeros((len(S),2))
    Srunning=0.0
    numpy.random.seed(seed=1234)
    nall=numpy.zeros((npix,len(S)))
    nallsum=numpy.zeros((2,len(S)))
    nallsum[0,:]=S
    SNP=numpy.zeros(npix*len(S))
    NP=numpy.zeros(npix*len(S))
    # For each pixel:
    for jpix in range(npix):
        np=numpy.zeros(len(S))
        # For each flux bin:
        for iflux in range(len(S)):
            # Poisson sample SKADS
            np[iflux]=numpy.random.poisson(lam=float(n[iflux]))
            #print n[iflux]
            #print np[iflux]
        # All the flux counts for this pixel
        nall[jpix,:]=np
        # Total flux for this pixel
        Sp[jpix]=(S*np).sum()
        # Separate pixel fluxes
        SNP[jpix*len(S):(jpix+1)*len(S)]=S*np
        #print Sp[jpix]
        # Keep running tally of highest-total-flux pixel
        if Sp[jpix] > Srunning:
            nprunning[:,0]=np
            nprunning[:,1]=S
            Srunning=Sp[jpix]
        #nallsum[1,]=np.sum()
        #print jpix,np.sum(),Sp[jpix]

    for iflux in range(len(S)):
        nallsum[1,iflux]=nall[:,iflux].sum()

    Sp *= Jy2muJy
    noise_vec=numpy.zeros(npix)
    #noise_vec=numpy.random.normal(0.0,SURVEY_NOISE,npix)
    Sn=numpy.zeros(npix)
    Sn=Sp+noise_vec

    SNP *= Jy2muJy

    #print Sp,Srunning
    #print nprunning
    #numpy.savetxt('snp.txt',SNP)
    #numpy.savetxt('sn.txt',S*N*Jy2muJy)
    return Sp,Srunning,nprunning,nallsum[1,:],Sn,SNP,S*N*Jy2muJy,S

#-------------------------------------------------------------------------------

@profile
def simtable(bins,a=None,Smin=None,Smax=None,N=None,A=None,seed=None,\
             noise=None,dump=False,output=None,verbose=False,version=2,\
             NLAWS_SIM=None,b=None,S0=None):

    """
    Ketron's Table 1 + simtable.pro

    This is working fine apart from the first and last bins
    (deriv. error?); actually I've now fixed that by manually fixing
    the endpoint derivatives to match IDL's DERIV.PRO

    noise is the sigma in uJy
    A (in sq. deg.) and the units are currently ignored
                        (as in Ketron's Table 1, despite what it says)
    Call as e.g.
        b=lumfunc.simtable(lumfunc.bins,a=-1.5,seed=1234,noise=10.0,dump=True)

    Note that N \propto C   ..   C \propto N (| alpha, Smin, Smax)

    """

    if N is None:
        if NLAWS_SIM==1:
            # Do this for AREA_SIM -> N total observed sources
            #N=ksRaw.sum() # N effectively sets C (or vice versa)
            # N or N2???? --- N I think, judging by Ketron's code
            N3=int(round(integrate.quad(lambda S:powerLawFuncS(\
                S,C_SIM,ALPHA_SIM,SMIN_SIM,SMAX_SIM,AREA_SIM/sqDeg2sr),\
                SMIN_SIM,SMAX_SIM)[0]))

            # Now just do this analytically:
            X1=SMIN_SIM; X2=SMAX_SIM
            #X1=110.3; X2=839.3
            A1=ALPHA_SIM+1
            N=int(round(C_SIM*((X2**A1-X1**A1)/A1) * AREA_SIM /sqDeg2sr))
            #N*=7.07
            N2=int(round(integrate.quad(lambda S:powerLawFuncErfsS(\
                S,NLAWS_SIM,C_SIM,ALPHA_SIM,-99.0,BETA_SIM,SMIN_SIM,SMAX_SIM,
                110.3,839.3,S0_SIM,GAMMA_SIM,S1_SIM,DELTA_SIM,S2_SIM,\
                NOISE_SIM,AREA_SIM/sqDeg2sr),\
                110.3,839.3)[0]))
            print N3,N,N2,C_SIM,ALPHA_SIM,SMIN_SIM,SMAX_SIM,NOISE_SIM,AREA_SIM
            #N=N2
        elif NLAWS_SIM==2:
            #N1=integrate.quad(lambda S:powerLawFuncS(S,C_SIM,ALPHA_SIM,SMIN_SIM,S0_SIM,AREA_SIM),SMIN_SIM,S0_SIM)[0]
            #N2=integrate.quad(lambda S:powerLawFuncS(S,D_SIM,BETA_SIM,S0_SIM,SMAX_SIM,AREA_SIM),S0_SIM,SMAX_SIM)[0]
            #N=N1+N2
            print C_SIM,ALPHA_SIM,D_SIM,BETA_SIM,SMIN_SIM,SMAX_SIM,S0_SIM,AREA_SIM
            N=int(round(integrate.quad(lambda S:powerLawFuncDoubleS(\
                S,C_SIM,ALPHA_SIM,D_SIM,BETA_SIM,SMIN_SIM,\
                SMAX_SIM,S0_SIM,AREA_SIM/sqDeg2sr),SMIN_SIM,SMAX_SIM)[0]))

    else:
        N4=C2N(C_SIM,ALPHA_SIM,SMIN_SIM,SMAX_SIM,AREA_SIM/sqDeg2sr)
        print "-> You've forced N: %i -> %i " % (N4,N)

    # It's N that's carried through


    # These should really be removed
    #if a is None: a=-1.5
    #if Smin is None: Smin=bins[0] # uJy
    #if Smax is None: Smax=bins[-1] # uJy

    #A1=ALPHA_SIM+1
    #CC=A1/(SMAX_SIM**A1-SMIN_SIM**A1)
    #N*= 1.0/CC
    #print N,CC,A1

    if NLAWS_SIM == 1:
        R = randomp(a,N,range_x=[X1,X2],seed=SEED_SIM,dump=None)
 
        #r = numpy.random.uniform(SMIN_SIM,SMAX_SIM,N)

        #fluxGen=buildPowerLawCDF(NLAWS_SIM,C_SIM/sqDeg2sr,ALPHA_SIM,-99.0,-99.0,\
        #                               SMIN_SIM,SMAX_SIM,-99.0,AREA_SIM,deltaS=0.01)
        #R = numpy.nan*numpy.ones(len(r))
        #for index,sample in enumerate(range(len(r))):
        #    R[index]=sampleBPL(fluxGen,r[index])

        #for index in range(len(r)):
        #    R[index]=powerLawFuncS(NLAWS_SIM,r[index]/1.0e6,C_SIM,ALPHA_SIM,
        #                           SMIN_SIM/1.0e6,SMAX_SIM/1.0e6,AREA_SIM)
        #    print r[index],R[index]

        #dn_by_ds=calculateDnByDs(bins,counts,verbose=verbose)
    elif NLAWS_SIM in [2,3,4]:
        N=70000
        fluxGen=buildPowerLawCDF(NLAWS_SIM,C_SIM,ALPHA_SIM,D_SIM,BETA_SIM,\
                                SMIN_SIM,SMAX_SIM,S0_SIM,GAMMA_SIM,\
                                S1_SIM,DELTA_SIM,S2_SIM,AREA_SIM,deltaS=0.01)

        if seed is not None:
            numpy.random.seed(seed=SEED_SIM)

        r = numpy.random.rand(N)

        # How to vectorize this loop?
        R=numpy.zeros(len(r))
        for index,sample in enumerate(range(len(r))):
            R[index]=sampleBPL(fluxGen,r[index])

    elif NLAWS_SIM == 0:
        N=72000
        skadsf='heywood/wilman_counts_1p4GHz_jz.txt'
        #R,w,x,y,Sn,SNP=sampleSKADS(skadsf=skadsf)
        R,w,x,NALLSUM,Sn,SNP,NORIG,S=sampleSKADS(skadsf=skadsf)
        #R=SNP
        #R=NORIG
        #R=NALLSUM
        #R=numpy.repeat(S*Jy2muJy,NALLSUM.astype('int64'),axis=0)
        R=Jy2muJy*10**numpy.genfromtxt('heywood/1sqdeg_0p02uJy.txt')
        numpy.ndarray.sort(R)
        Rmin=SMIN_SKADS; Rmax=SMAX_SKADS
        iRmin,Rmin=find_nearest(R,Rmin)
        iRmax,Rmax=find_nearest(R,Rmax)
        print iRmin,Rmin,iRmax,Rmax
        R=R[iRmin:iRmax]
        N=len(R)
        if NSKADS is not None:
            R=numpy.random.choice(R,size=NSKADS,replace=False)
            N=len(R)
        print 'NSKADS = %i' % N
        #print R.max(),R.min()

    elif NLAWS_SIM == -1: # Not used

        # Build the CDF
        skadsf='heywood/wilman_counts_1p4GHz_jz.txt'
        skads=numpy.genfromtxt(skadsf)
        x=skads[:,0] # S
        y=numpy.power(skads[:,0],-2.5)*skads[:,3] # dnds
        z=skads[:,2] # N
        numpy.savetxt('y.txt',y)

        #lookup = interp1d(x,y)
        #values=numpy.array([lookup(S) for S in x])

        xmin=x[0]; xmax=x[-1]
        xmin=0.01e-6; xmax=85.0e-6
        ixmin,xmin=find_nearest(x,xmin)
        ixmax,xmax=find_nearest(x,xmax)
        print xmin,xmax
        print len(x[ixmin:ixmax])
        #Nskads_sqdeg_cutdown=skads[ixmin:ixmax,2].sum()
        #vc=y.cumsum()/y.cumsum()[-1]
        Nskads_sqdeg=z.sum()
        print  Nskads_sqdeg
        Nskads_sqdeg *= sqDeg2sr
        N=int(Nskads_sqdeg)
        N_in_range=z[ixmin:ixmax].sum()
        print z[ixmin:ixmax].cumsum()
        N_in_range *= sqDeg2sr
        N*=10.0
        print 'N in SKADS sim',N,N_in_range
        N=N
        #N=int(N_in_range)
        #x=x[ixmin:ixmax]; z=z[ixmin:ixmax]; y=y[ixmin:ixmax]
        zz=z.cumsum()/z.cumsum()[-1]
        yy=y.cumsum()/y.cumsum()[-1] # i.e. =cumsum/sum
        print zz.shape,z.sum(),z.sum()
        #print z.cumsum()[-1],z.sum()
        #print len(vc)
        #numpy.savetxt('vc.txt',vc)
        #numpy.savetxt('zz.txt',zz)
        #numpy.savetxt('yy.txt',yy)
        #print min(vc),max(vc)
        print xmin,xmax
        fluxGen=interp1d(zz,x)

        if seed is not None:
            numpy.random.seed(seed=SEED_SIM)
        r = numpy.random.rand(N)
        print r.min(),r.max()
        R=numpy.zeros(len(r))
        #for i in numpy.linspace(0.1,1.0,110)[:-1]:
        #    print Jy2muJy*sampleBPL(fluxGen,i)

        # How to vectorize this loop (quickly)?
        for index,sample in enumerate(range(len(r))): # overkill!
            R[index]=Jy2muJy*sampleBPL(fluxGen,r[index])
            if R[index] < 0.05 or R[index] > 80.0:
                continue
            #R[index]='NaN'
        #R = R[~numpy.isnan(R)]
        print len(R)

    # Dump noiseless fluxes to file
    if dump is not None:
        dumpf=dump
        numpy.savetxt(dumpf,R)
        print 'Draws (noiseless) are in %s' % dumpf
    # Add noise
    if SIM_DO_CAT_NOISE: # switch off noise temporarily if want
        numpy.random.seed(seed=SEED_SIM)
        R+=numpy.random.normal(0.0,noise,N)

    # Dump noisy fluxes to file
    if dump is not None:
        noisydumpf='%s_noisy.txt' % dumpf.split('.')[0]
        numpy.savetxt(noisydumpf,R)
        print 'Draws (noisy) are in %s' % noisydumpf
        print 'Minimum flux in catalogue = %f' % R.min()

    # Bin up the fluxes
    counts=numpy.histogram(R,bins=bins)[0]

    print '-> %i/%i objects observed in total (bin scattering)\n' % (counts.sum(),N)

    # Now calculate dn/ds, dn/ds_eucl
    # bins assumed to be in uJy
    # counts are per unit SURVEY_AREA (AREA_SIM)
    # dn_by_ds(_eucl) are in standard units
    idl_style=False
    dn_by_ds=calculateDnByDs(1.0e-6*bins,counts,eucl=False,idl_style=idl_style)
    dn_by_ds_eucl=calculateDnByDs(1.0e-6*bins,counts,eucl=True,idl_style=idl_style)
    dn_by_ds_errs=calculateDnByDs(1.0e-6*bins,counts,errors=True,idl_style=idl_style)
    #print dn_by_ds
    #print dn_by_ds_eucl

    median_bins=medianArray(bins) # uJy
    if output is not None:
        outputf=output
        #data=numpy.zeros((len(bins),3))
        s=open(outputf,'w')
        if version < 2:
            header='# bin_median ksRaw ksNoisy'
        else:
            #header='# bin_low_uJy bin_high_uJy bin_median_uJy Ns_tot_obs dnds_ dnds_eucl delta_dnds_eucl_lower delta_dnds_eucl_upper corr Ngts dNgts_lower dNgts_upper'
            header='# bin_low_uJy bin_high_uJy bin_median_uJy Ns_tot_obs dnds_srm1Jym1 dnds_eucl_srm1Jy1p5 delta_dnds_eucl_lower_srm1Jy1p5 delta_dnds_eucl_upper_srm1Jy1p5 corr Ncgts_degm2 dNcgts_lower_degm2 dNcgts_upper_degm2'
        s.write('%s\n'%header)
        print header
        for ibin in range(nbins-1):
            #data[ibin,:]=(bins[ibin],bins[ibin+1],dn[ibin])
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
                   round(sqrt(counts[ibin:].sum()*1.00/AREA_SIM)),\
                   round(sqrt(counts[ibin:].sum()*1.00/AREA_SIM)))
            s.write('%s\n'%line)
            print line
        print counts.sum()
        s.close()


    return counts

#-------------------------------------------------------------------------------

# Vestigial part of simtable function

    #A=10.0 # sq. deg.
    #N=ksRaw.sum()
    #print N*A
###    N *= A/REF_AREA # / sq. deg.
    # This is NOT how to do it
    #DeltaS=Smax-Smin # uJy
    #Scentre=(Smin+Smax)/2.0 # uJy
    #dn_over_ds=N/DeltaS # / uJy / sq. deg.


    #a=-1.5
    # Generate the random draws
    R = randomp(a,N,range_x=[Smin,Smax],seed=seed,dump=dump)
    #R = powerlaw.rvs(a, size=100)
    bRaw=numpy.histogram(R,bins=bins)
    dn_raw=bRaw[0]
    dn_by_ds_raw=calculateDnByDs(bins,bRaw[0],verbose=verbose)
    log10_dn_by_ds_raw = numpy.log10(dn_by_ds_raw)


    # Now add noise and bin again
    if noise is not None:
        print '(Please wait a moment...)'
        #for isamp in range(len(R)):
        #    R[isamp] = numpy.random.normal(R[isamp],noise,1)
        # There's probably a better way to write this:
        R=map(lambda y:(numpy.random.normal(y,noise,1)),R)

    b=numpy.histogram(R,bins=bins)

    print 'There are %i/%i objects in total (bin scattering)\n' % (b[0].sum(),N)

    #print 'b',b
    dn_by_ds=calculateDnByDs(bins,b[0],verbose=verbose)
    # ^^^^^^^^^^ now in ^^^^^^^^^^^^
    dn=b[0]
    #ds=numpy.absolute(numpy.gradient(median_bins))
    #print 'dn',dn
    #print 'ds',ds
    #dn_by_ds = dn/ds
    log10_dn_by_ds = numpy.log10(dn_by_ds)

    median_bins=medianArray(bins)

    if verbose:
        print 'log10(dn/ds)',log10_dn_by_ds
        print '# bin_low bin_high bin_median dn_raw log10_dn_by_ds_raw dn_noisy log10_dn_by_ds_noisy'
        for ibin in range(nbins-1):
            print '%6.4f %6.4f %6.4f %i %6.4f %i %6.4f' % \
              (bins[ibin],bins[ibin+1],median_bins[ibin],
               dnRaw[ibin],log10_dn_by_ds_raw[ibin],dn[ibin],log10_dn_by_ds[ibin])

    if output is not None:
        outputf=output
        #data=numpy.zeros((len(bins),3))
        s=open(outputf,'w')
        header='# bin_median ksRaw ksNoisy'
        for ibin in range(nbins-1):
            #data[ibin,:]=(bins[ibin],bins[ibin+1],dn[ibin])
            line='%f %i %i\n' % (median_bins[ibin],dn_raw[ibin],dn[ibin])
            s.write(line)
        s.close()

    return dn_by_ds

#-------------------------------------------------------------------------------

@profile
def buildCDF(func,S1,S2,deltaS):
    """
    Following Russell's prescription [email of 24.1.14]
    Return an array useable for a mock generation for an arbitrary function
    """

    x=numpy.arange(S1,S2+deltaS,deltaS)
    y = numpy.array([func(ix) for ix in x])
    f = interp1d(x,y)

    return f

#-------------------------------------------------------------------------------

@profile
def buildPowerLawCDF(nlaws,C,alpha,D,beta,Smin,Smax,S0,gamma,S1,delta,S2,area,deltaS=0.01):

    """
    I think there may be an added layer of complication used here
    But this only has to be called once for the simulation
    So leave it for now...
    """

    lookup=buildCDF(lambda S:powerLawFuncWrap(\
                nlaws,S,C,alpha,-99.0,beta,Smin,Smax,S0,gamma,S1,delta,S2,area),Smin,Smax,deltaS)
                
    Ss=numpy.arange(Smin,Smax+deltaS,deltaS)
    values=numpy.array([lookup(S) for S in Ss])
    vc=values.cumsum()/values.cumsum()[-1]

    return interp1d(vc,Ss)

#-------------------------------------------------------------------------------

@profile
def sampleBPL(lookup,x):
    """
    For a given sample U[0,1], return a sample flux from a (possibly broken)
    power law built using buildPowerLawCDF
    """
    
    return float(lookup(x))
    
#-------------------------------------------------------------------------------

from scipy.interpolate import UnivariateSpline

@profile
def poles(nlaws,Smin,Smax,noise=None):

#def calculateI3(C,alpha,Smin,Smax,area,noise=None,dump=None,\
#                verbose=False,nlaws=1,D=None,beta=None,S0=None):
    return



#-------------------------------------------------------------------------------

@profile
def counts2TemperatureArray(synthBeamSolidAngleSr,fluxes,counts,Smin,Smax,freq=None):
    """
    Integrate an array of counts, e.g. Vernstrom
    """

    #mask=numpy.where(fluxes<=Smax) and numpy.where(fluxes>=Smin)
    mask=numpy.where(fluxes<=Smax)
    print Smin,Smax,fluxes[mask][0],fluxes[mask][-1]
    Tb=1.0e26*synthBeamSolidAngleSr\
      *integrate.trapz(fluxes[mask]*counts[mask],fluxes[mask])

    Tb*=(2.0*kb*freq**2)/clight**2

    return Tb

#-------------------------------------------------------------------------------

@profile
def counts2Temperature(nlaws,C,alpha,beta,Smin,Smax,S0,gamma,S1,delta,S2,freq=None):

    """
    Vernstrom 2014 equation (14)
    For S in Jy,
    Returns Tb in K
    """
    DUMMY=1.0
    Tb=1.0e26*radioSynthOmegaSr*integrate.quad(lambda S:S*powerLawFuncWrap(nlaws,S,\
        C,alpha,-99.0,beta,Smin,Smax,S0,gamma,S1,delta,S2,DUMMY),Smin,Smax)[0]
    Tb*=(2.0*kb*freq**2)/clight**2

    return Tb

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
        #print S,S0,S1,alpha,beta,gamma
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
        #print S,S0,S1,S2,alpha,beta,gamma,delta
        n = C  * S0**(alpha-beta) * S1**(beta-gamma) * \
                 S2**(gamma-delta)* (S**delta) * area

    return n

#-------------------------------------------------------------------------------

@profile
def piecewise5(x,n0,posns=[1.0,2.0,3.0],alpha=[-2.0,-1.5]):
    """
    Not working yet - produces jumps!
    """

    assert(len(alpha)==len(posns)-1)
    if not (posns[0] < float(x) <= posns[-1]): return 0.0

    na=len(alpha)
    ithis = next(i for i in xrange(na) if (posns[i] < float(x) <= posns[i+1]))
    if float(x) in posns:
        print x-1,x,x+1,ithis,posns[ithis],posns[ithis+1],[i for i in xrange(ithis+1)]
        print n0*numpy.cumprod([(float(x)**alpha[i]*posns[i]**(alpha[i]-alpha[i+1])) for i in xrange(ithis)])
        #print ithis,posns[ithis],float(x),posns[ithis+1],[i for i in xrange(ithis+1)]
    if ithis==0: return n0*float(x)**alpha[0] # annoying hack
    n=n0*numpy.cumprod([(float(x)**alpha[i]*posns[i]**(alpha[i]-alpha[i+1])) for i in xrange(ithis)]).tolist()[-1]
    return n

#-------------------------------------------------------------------------------

@profile
def piecewise4(x,n0,posns=[1.0,2.0,3.0],alpha=[-2.0,-1.5]):
    """
    This is a linearized version of
    n = C . S ^ alpha
    n *= ......
    giving a linear piecewise power-law s.t.

      { C . S ^ alpha[0], posns[0] < n < posns[1]
    n={ C . S ^ ....    , posns[1] < n < posns[2]
      { C . S ^ ....    , posns[2] < n < posns[...]
      { 0                 otherwise

    3.8.14 [see notebook and stackoverflow]
    """

    assert(len(alpha)==len(posns)-1)

    if not (posns[0] < x <= posns[-1]): return 0.0

    ln_n = log(n0) + alpha[0]*log(float(x))

    np=len(posns)

    for ip in range(np):
        if posns[ip] < x <= posns[ip+1]: return exp(ln_n)
        ln_n += (alpha[ip+1]-alpha[ip])*(log(float(x))-log(posns[ip+1]))
    return exp(ln_n)

#-------------------------------------------------------------------------------

@profile
def piecewise3(x,n0,posns=[1.0,2.0,3.0],alpha=[-2.0,-1.5]):

    """
    """

    assert(len(alpha)==len(posns)-1)

    if x <= posns[0] or x > posns[-1]: return 0.0

    n=n0*x**alpha[0]

    np=len(posns)
    #for ip in range(np):
    #    if posns[ip] < x <= posns[ip+1]: return n
    #    n *= (posns[ip+1]/float(x))**(alpha[ip]-alpha[ip+1])
    #return n

    ns=[(posns[ip+1]/float(x))**(alpha[ip]-alpha[ip+1]) for ip in range(np-2)]
    print 'unfinished'
    return ns[0]


#-------------------------------------------------------------------------------

@profile
def piecewise2(x,n0,posns=[1.0,2.0,3.0],alpha=[-2.0,-1,5]):

    if x <= posns[0] or x > posns[-1]: return 0.0

    n=n0*x**alpha[0]
    np=len(posns)
    for ip in range(np):
        if posns[ip] < x <= posns[ip+1]: return n
        n *= (posns[ip+1]/float(x))**(alpha[ip]-alpha[ip+1])
    return n

#-------------------------------------------------------------------------------

@profile
def piecewise(x,ampl,posns=[1.0,2.0,3.0],alpha=[-2.0,-1,5]):

    if x <= posns[0] or x > posns[-1]: return 0.0

    n=ampl*x**alpha[0]

    np=len(posns)
    for ip in range(np):
        print posns[ip],x,posns[ip+1],n
        if posns[ip] < x <= posns[ip+1]: return n
        n *= x**(alpha[ip]-alpha[ip+1])
    return n

#-------------------------------------------------------------------------------

@profile
def powerLawFuncArb(S,nlaws,params,area):

    """
    All fluxes in Jy by convention
    
    Order of parameters
    0 C
    1 alpha
    2 Smin
    3 Smax
    4 beta
    5 S0
    6 gamma
    7 S1
    8 delta
    9 S2
    10 ...
    11 ...
    """

    C=alpha=Smin=Smax=beta=S0=gamma=S1=delta=S2=-99.0
    if nlaws==1:
        (C,alpha,Smin,Smax)=params[:4]
        #print C,alpha,Smin,Smax
        return powerLawFuncS(S,C,alpha,Smin,Smax,area)
    elif nlaws==2:
        (C,alpha,Smin,Smax,beta,S0)=params[:6]
        return powerLawFuncDoubleS(S,C,alpha,-99.0,beta,Smin,Smax,S0,area)
    elif nlaws==3:
        (C,alpha,Smin,Smax,beta,S0,gamma,S1)=params[:8]
        return powerLawFuncTripleS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,area)
    elif nlaws==4:
        (C,alpha,Smin,Smax,beta,S0,gamma,S1,delta,S2)=params[:10]
        return powerLawFuncQuadS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,delta,S2,area)

    return

#-------------------------------------------------------------------------------

@profile
def computeSmin(nlaws,C,alpha,beta,Smax,S0,gamma,S1,delta,S2,area,NTHRESH):

    #func = powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Smin,Smax,S0,gamma,S1,delta,S2,area)
    Smin=0.01
    yS=numpy.arange(1.0,S1,0.1)[::-1] # reverse
    #print yS
    N_S1_SMAX=integrate.quad(lambda S:powerLawFuncWrap(nlaws,S,C,alpha,-99.0,beta,Smin,Smax,S0,gamma,S1,delta,S2,area),S1,Smax)[0]
    N_Si_SMAX=N_S1_SMAX
    for i,Si in enumerate(yS):
        N_Si_SMAX+=integrate.quad(lambda S:powerLawFuncWrap(nlaws,S,C,alpha,-99.0,beta,Smin,Smax,S0,gamma,S1,delta,S2,area),Si,yS[i+1])[0]
        print i,Si,yS[i+1],integrate.quad(lambda S:powerLawFuncWrap(nlaws,S,C,alpha,-99.0,beta,Smin,Smax,S0,gamma,S1,delta,S2,area),yS[i+1],Si)[0]
        if N_Si_SMAX >= NTHRESH:
            return Si

    return 0

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

    # This was my original double power law
    if S <= Smin or S > Smax:
    #if False:
        return 0.0
    elif Smin < S <= S0:
    #elif S <= S0:
        n = C * (S ** alpha) * area
    elif S0 < S <= Smax:
    #elif S > S0:
        n = D * (S ** beta) * area
    return n

#-------------------------------------------------------------------------------

@profile
def powerLawFuncS(S,C,alpha,Smin,Smax,area):
    """
    Ketron's equation (9)
    """

    if dnds0 and (S < Smin or S > Smax):
        return 0.0
    else:
        n = C * (S ** alpha) * area
        #print n, C, S, alpha, area
        return n

#-------------------------------------------------------------------------------

@profile
def countFunction(model):

    """
    Generalized version of Ketron's equation (5) + (9), for:
    power laws
    modified power laws
    polynomials
    nodal models
    """

    print model.name

    return realisation

#-------------------------------------------------------------------------------

@profile
def powerLawFuncArbErfsS(S,nlaws,params,Sbinlow,Sbinhigh,sigma,area):
    """
    All fluxes in Jy
    sigma in Jy
    """

    Smin=params[2]; Smax=params[3]
    if S < Smin or S > Smax: return 0.0
    assert(sigma>1.0e-50), 'Noisy data now required %f' % sigma

    erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*sigma)) - erf((S-Sbinhigh)/(sqrtTwo*sigma)))
    #print 'a',S,params,Sbinlow,Sbinhigh,sigma,erfs
    return erfs * powerLawFuncArb(S,nlaws,params,area)


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

    if ssigma > 1.0e-50:
        erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))
        #erfs *= 0.5
        # Should the next line be min/max or binlow/binhigh???
        # min/max - parameter inferences way off otherwise
        #return erfs * powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Sbinlow,Sbinhigh,S0,area)
        return erfs * powerLawFuncWrap(nlaws,S,C,alpha,D,beta,\
                                       Smin,Smax,S0,gamma,S1,delta,S2,area)
    else: # Noise-free(!!)
        # Should the next line be min/max or binlow/binhigh???
        # THIS NEEDS TO BE THOROUGHLY CHECKED
        return powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Sbinlow,Sbinhigh,\
                         S0,gamma,S1,delta,S2,area)

#-------------------------------------------------------------------------------


x = numpy.linspace(-5.0,5.0,100)
y=erf(x)
f = interp1d(x, y)
def timeErf2(xx):
    return f(xx)

def timeErf1(xx):
    return erf(xx)

#-------------------------------------------------------------------------------



pri=Priors()

@profile
def myprior(cube,ndim,nparams):
    """
    ##model 0: vary no parameters
      model 1: vary all 4 [or 6] parameters
    ##model 2: vary C and alpha only
    ##model 3: vary Smin and Smax only
    ##model 4: vary alpha, Smin and Smax only
    ##model 10: vary C, alpha, beta, Smin, Smax, S0
    
    Try to get scipy.stats intrinsics working (faster? but disagrees!?)
    #PriorU0 = stats.uniform(1.0,100.0)
    #PriorU1 = stats.uniform(-2.5,-0.1)
    #PriorU2 = stats.uniform(0.01,5.00)
    #PriorU3 = stats.uniform(1.0*sigma,5.0*sigma)
    cube[0]=PriorU0.ppf(cube[0])
    cube[1]=PriorU1.ppf(cube[1])
    cube[2]=PriorU2.ppf(cube[2])
    cube[3]=PriorU3.ppf(cube[3])
    """
    
    # Vary no parameters
    #if model == 0:
    #    cube[0]=pri.DeltaFunctionPrior(cube[0],40.0,40.0)
    #    cube[1]=pri.DeltaFunctionPrior(cube[1],-2.5,-2.5)
    #    cube[2]=pri.DeltaFunctionPrior(cube[2],1.0,1.0)
    #    cube[3]=pri.DeltaFunctionPrior(cube[3],20.0,20.0)

    # Vary all 4 parameters - C, alpha, Smin, Smax [, beta, S0]
    #elif model == 1:
    cube[0]=pri.GeneralPrior(cube[0],C_PRIOR,C_MIN,C_MAX) # C / 
        #if C_PRIOR=='U':
        #    cube[0]=pri.UniformPrior(cube[0],C_MIN,C_MAX) # C / /Jy /sr
        #elif C_PRIOR=='LOG':
        #    cube[0]=pri.LogPrior(cube[0],C_MIN,C_MAX)
        #cube[0]=40.0
    cube[1]=pri.UniformPrior(cube[1],ALPHA_MIN,ALPHA_MAX) # alpha
        #cube[1]=-1.5
        #cube[2]=pri.UniformPrior(cube[2],0.01e-6,5.00e-6) # Smin / Jy
        #cube[2]=pri.UniformPrior(cube[2],0.01,5.00) # Smin / uJy
        #cube[3]=pri.UniformPrior(cube[3],1.0*sigma,5.0*sigma) # Smax / uJy
    cube[2]=pri.UniformPrior(cube[2],SMIN_MIN,SMIN_MAX) # Smin / uJy
    cube[3]=pri.UniformPrior(cube[3],SMAX_MIN,SMAX_MAX) # Smax / uJy
    if nlaws > 1:
        #cube[4]=pri.UniformPrior(cube[4],D_MIN,D_MAX) # D / /Jy /deg^2
        cube[4]=pri.UniformPrior(cube[4],BETA_MIN,BETA_MAX) # beta
        cube[5]=pri.UniformPrior(cube[5],S0_MIN,S0_MAX) # S0
    if nlaws > 2:
        cube[6]=pri.UniformPrior(cube[6],GAMMA_MIN,GAMMA_MAX) # beta
        cube[7]=pri.UniformPrior(cube[7],S1_MIN,S1_MAX) # S0
    if nlaws > 3:
        cube[8]=pri.UniformPrior(cube[8],DELTA_MIN,DELTA_MAX) # beta
        cube[9]=pri.UniformPrior(cube[9],S2_MIN,S2_MAX) # S0

    if floatNoise:
        cube[nparams-1]=pri.UniformPrior(cube[nparams-1],NOISE_MIN,NOISE_MAX) # sigma

    # Fix Smin and Smax - vary just C and alpha
    #elif model == 2:
    #    cube[0]=pri.UniformPrior(cube[0],1.0,100.0) # C / /Jy /deg^2
    #    cube[1]=pri.UniformPrior(cube[1],-2.5,-0.1) # alpha
    #    cube[2]=pri.DeltaFunctionPrior(cube[2],1.0,1.0)
    #    cube[3]=pri.DeltaFunctionPrior(cube[3],20.0,20.0)

    # Fix C and alpha - vary just Smin and Smax
    #elif model == 3:
    #    cube[0]=pri.DeltaFunctionPrior(cube[0],40.0,40.0)
    #    cube[1]=pri.DeltaFunctionPrior(cube[1],-1.5,-1.5)        
    #    cube[2]=pri.UniformPrior(cube[2],0.01,2.00) # Smin / uJy
    #    cube[3]=pri.UniformPrior(cube[3],10.0,30.0) # Smax / uJy

    # Fix C - vary just alpha, Smin and Smax
    #elif model == 4:
    #    cube[0]=pri.DeltaFunctionPrior(cube[0],40.0,40.0)
    #    cube[1]=pri.UniformPrior(cube[1],-2.5,-0.1)        
    #    cube[2]=pri.UniformPrior(cube[2],0.01,2.00) # Smin / uJy
    #    cube[3]=pri.UniformPrior(cube[3],10.0,30.0) # Smax / uJy

    # Fix alpha - vary just C, Smin and Smax
    #elif model == 5:
    #    cube[0]=pri.UniformPrior(cube[0],1.0,100.0) # C / /Jy /deg^2
    #    cube[1]=pri.DeltaFunctionPrior(cube[1],-1.5,-1.5)        
    #    cube[2]=pri.UniformPrior(cube[2],0.01,2.00) # Smin / uJy
    #    cube[3]=pri.UniformPrior(cube[3],10.0,30.0) # Smax / uJy

    # Fix Smin - vary just C, alpha and Smax
    #elif model == 6:
    #    cube[0]=pri.UniformPrior(cube[0],1.0,100.0) # C / /Jy /deg^2
    #    cube[1]=pri.UniformPrior(cube[1],-2.5,-0.1)        
    #    cube[2]=pri.DeltaFunctionPrior(cube[2],1.00,1.00) # Smin / uJy
    #    cube[3]=pri.UniformPrior(cube[3],10.0,30.0) # Smax / uJy

    # Fix Smax - vary just C, alpha and Smin
    #elif model == 7:
    #    cube[0]=pri.UniformPrior(cube[0],1.0,100.0) # C / /Jy /deg^2
    #    cube[1]=pri.UniformPrior(cube[1],-2.5,-0.1)        
    #    cube[2]=pri.UniformPrior(cube[2],0.01,2.00) # Smin / uJy
    #    cube[3]=pri.DeltaFunctionPrior(cube[3],20.0,20.0) # Smax / uJy

    # Double power law
    #elif model == 10:
    #    cube[0]=pri.UniformPrior(cube[0],1.0,100.0) # C / /Jy /deg^2
    #    cube[1]=pri.UniformPrior(cube[1],-2.5,-0.1)        
    #    cube[2]=pri.UniformPrior(cube[2],0.01,2.00) # Smin / uJy
    #    cube[3]=pri.UniformPrior(cube[3],10.0,30.0) # Smax / uJy
    #    #cube[4]=pri.UniformPrior(cube[4],1.0,100.0) # D / /Jy /deg^2
    #    cube[4]=pri.UniformPrior(cube[4],-2.5,-0.1)
    #    cube[5]=pri.UniformPrior(cube[5],0.01,30.0) # S0 / uJy

    if loud:
        if True:
        #    # http://stackoverflow.com/questions/16578652/threading-timer
            updatePlot(plothandle,cube)
        else:
            pass

    return

#-------------------------------------------------------------------------------

#def calculateIDouble(C,alpha,Smin,Smax,Sbinlow,Sbinhigh,sigma):
#    I=integrate.dblquad(powerLawFuncErfsS(S,C,alpha,Smin,Smax,\
#                        Sbinlow,Sbinhigh,sigma),\
#                        Smin,Smax,lambda S:Sbinlow,lambda S:Sbinhigh,args=())[0]
#    return I

#-------------------------------------------------------------------------------

@profile
def calculateI3(C,alpha,Smin,Smax,area,noise=None,dump=None,\
                verbose=False,nlaws=1,D=None,beta=None,S0=None,\
                gamma=-99.0,S1=-99.0,delta=-99.0,S2=-99.0):
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


    #II = numpy.zeros(nbins-1)
    ## Set up the numerical derivative
    #resolution=1000
    #ds = (Smax - Smin) / resolution
    #for ibin in xrange(nbins-1):
    ## This is the integral from Smin to Smax
    ## Step from Smin -> Smax
    #    I=0.0
    #    for istep in xrange(resolution):
    #        Si = Smin + (istep + 0.5) * ds # + -> - because of counting from 0 not 1
    #        erfs = erf((Si-bins[ibin])/(sqrt(2.0)*noise)) - erf((Si-bins[ibin+1])/(sqrt(2.0)*noise))
    #        I += area*C*Si**alpha * erfs * ds
    #    II[ibin]=I

    # I checked this is the same as the numerical integration above
    # (to 1 dec. place.)
    II = numpy.zeros(nbins-1)
    II2 = numpy.zeros(nbins-1)
    # Switch to quadrature integration
    # This is the integral from Smin to Smax
    # Put C scalings here:
    #C *= 10**(1.0/(6.0*alpha))
    #C *= 10**(6.0*alpha)
    #C *= 10**6.0 # Get junk out unless C is of this order (uJy v. Jy)
    for ibin in xrange(nbins-1):
        # This is switched off in order to fix the units of C, but for
        # a long time this was the version in use (1401NNx)
        #II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,area),Smin,Smax)[0]*1000.0

        # I'm trying to investigate the units of the power law in
        # case it affects C
        # **** !!!! This is the wrong place to do it --- the loop compounds!!!!
        # Tried adding a factor of SURVEY_AREA here:
        #print C,alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,area
        #raw_input('pause')
        # 1.0 was = area but I've hard-coded this now
        #C *= 1.0/sqDeg2sr
        #C *= 1.0e6 <- seems to work ok at 140207, but error bars too small
        #C *= 10**(1.0/(6.0*alpha))
        #**** All the quantities are in uJy, except C is in uJy^-1 sr^-1
        #**** alpha is dimensionless
        #**** II needs to be in same area units as ksNoisy (/SURVEY_AREA now)
        sqDeg2srr=sqDeg2sr
        #sqDeg2srr=1.0
        if nlaws == 1:
            # Could 'lambda: S' be deleted - speed-up..?
            # Was Smin,Smax -> Sbinlow,Sbinhigh
            # Need to sort out this maths C expression - II or II2??
            #II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C*10**(1.0/(6.0*alpha)),alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,1.0),bins[ibin],bins[ibin+1])[0]
            # Checked II and II2 are now the same:
            II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C*1.0e-6*10**(-6.0*alpha),alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,sqDeg2srr*SURVEY_AREA),Smin,Smax)[0]
            II2[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2srr*SURVEY_AREA),Smin/1.0e6,Smax/1.0e6)[0]
            #print bins[ibin],bins[ibin+1],Smin,Smax,II[ibin],II2[ibin],II[ibin]/II2[ibin]

            #print C,II[ibin],' '
            #II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,1.0),Smin,Smax)[0]
            #print II[ibin]

        elif nlaws in [2,3,4]:
            # XXXXX NEED TO FIX INTEGRATION LIMITS XXXXX -> Sbinlow,Sbinhigh
            II[ibin]=integrate.quad(lambda S:powerLawFuncErfsS(S,nlaws,C,alpha,D,beta,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,S0/1.0e6,gamma,S1/1.0e6,delta,S2/1.0e6,noise/1.0e6,sqDeg2srr*SURVEY_AREA),Smin/1.0e6,Smax/1.0e6)[0]

            #print ibin,bins[ibin],bins[ibin+1],Smin,Smax,II[ibin],II2[ibin],II[ibin]/II2[ibin]
        #print II.sum()*1000.0
    #print C,alpha,Smin,Smax
    #sys.exit(0)

    return II

#-------------------------------------------------------------------------------

@profile
def realiseData(nlaws,params,area,noise,bins):
    """
    Rewrite of calculateI3
    Bins in uJy
    sigma in uJy
    parameter fluxes in Jy
    """
    print params,noise,bins[0]
    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        II[ibin]=integrate.quad(lambda S:powerLawFuncArbErfsS(S,nlaws,params,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2sr*SURVEY_AREA),bins[ibin]/1.0e6,bins[ibin+1]/1.0e6)[0]
    print II
    return II

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

    dnds00=True
    if dnds00 and (Si < Smin or Si > Smax):
        #print '***Si',Smin,Smax
        #if False:
        return 0.0

    erfs = erf((Si-Sbinlow)/(sqrt(2.0)*noise)) - erf((Si-Sbinhigh)/(sqrt(2.0)*noise))
    erfs *= 0.5
    n = C * Si**alpha * erfs * area

    return n

#-------------------------------------------------------------------------------

@profile
def calculateISim(C,alpha,Smin,Smax,area,noise=None,seed=None,dump=None,\
                verbose=False):
    """
    Generate a binned array of noisy (Po + G) source counts

    lumfunc.simtable now makes use of this

    N total is set by the model parameters
    N is then used in the flux draws
    Then noise is added
    And the fluxes are binned up
    Hey presto!

    Was called calculateI3

    """

    # Initialize the sampling seed, if requested
    if seed is not None: numpy.random.seed(seed=SEED_SIM)

    # [get_nsource.pro; qsimp on dnds_model]
    # C needs to be multiplied up here (do this a few lines below instead?)

    # QSIMP, 'dnds_model', SMIN, SMAX, N

    # I've checked that (modulo a factor of 1000) this is 621114 for
    # 40.0,-1.5,1.0,20.0 [get_nsource + qsimp]
    # The factor of 1000 is related to the area (N \propto C * area) as far as I can tell
    # This is just C -> N, which we need for the generators below:
    N=integrate.quad(lambda S:powerLawFuncS(S,C,alpha,Smin,Smax,area),\
                     Smin,Smax)[0]*1000.0

    # This next line should mirror Ketron's randomn(seed,1,POISSON=trueN)
    # JZ NEED TO CHECK [randomn]
    Nold=N
    N = numpy.random.poisson(Nold)
    #print Nold,N

    # Generate the random draws (power-law distribution) [randomp]
    dump=None
    R = randomp(alpha,N,range_x=[Smin,Smax],seed=seed,dump=dump)
    #R = powerlaw.rvs(a, size=100)
    #print 'R', len(R)
    #numpy.savetxt('Rraw.txt',R)

    # Now add noise and - only then - bin
    if noise is not None:
        #print '(Please wait a moment...)'
        #for isamp in range(len(R)):
        #    R[isamp] = numpy.random.normal(R[isamp],noise,1)
        R+=numpy.random.normal(0.0,noise,N)
        #R=map(lambda y:(numpy.random.normal(y,noise,1)),R)

        #numpy.savetxt('Rnoisy.txt',R)

    #print 'R', len(R)

    # Would be good to write R out here sometimes

    #print bins
    # Now bin the noisy fluxes (so we should lose some sources
    # < Smin and > Smax)
    counts=numpy.histogram(R,bins=bins)[0]

    #There are about 336000 sources, so the dn is roughly spot on for
    # 40.0,-1.5,1.0,20.0
    #print 'c', counts, counts.sum()

    # Now differentiate to get dN/dS for each bin
    # bins, counts -> dN/dS
    dn_by_ds_realn=calculateDnByDs(bins,counts,verbose=verbose)
    #log10_dn_by_ds_realn = numpy.log10(dn_by_ds_realn)
    
    #print 'There are %i/%i objects in total (bin scattering)\n' % (brealn[0].sum(),N)

    #print 'b',b
    #print 'b',bins

    #print int(dn_by_ds_realn.sum())
    #print dn_by_ds_realn
    return dn_by_ds_realn

    # Now we need to take dN/dS for this realisation and integrate it
    # (eqn 5) between Smin and Smax, taking account of the erfs factor

    # We do this for all bins simultaneously
    II = numpy.zeros(nbins-1)

    # Set up the numerical derivative
    resolution=1000
    ds = (Smax - Smin) / resolution
    for ibin in xrange(nbins-1):
    # This is the integral from Smin to Smax
    # Step from Smin -> Smax
        I=0.0
        for istep in xrange(resolution):
            Si = Smin + (istep + 0.5) * ds # + -> - because of counting from 0 not 1
            erfs = 0.5*(erf((Si-bins[ibin])/(sqrt(2.0)*noise)) - erf((Si-bins[ibin+1])/(sqrt(2.0)*noise)))
            I += dn_by_ds_realn[ibin] * erfs * ds
        II[ibin]=I

    #There are about 373500 sources, so I *think* this is working right
    # 40.0,-1.5,1.0,20.0
    #print II.sum()

    # See if this constraint helps...
    #II *= N_SIM / N
    #print C,N_SIM,II.sum(),Nold,N
    # Next line is the most successful yet....
    #II *= N/II.sum()
    #II *= 775527.0/II.sum()

    return II


#-------------------------------------------------------------------------------

def testIntegral(x1,x2,resn):
    """
    """
    JJ=0.0
    dx = (x2-x1)/resn
    for istep in xrange(resn):
        xi=x1 + (istep-0.5)*dx
        JJ += dx*(xi**3)/3.0
    return JJ

@profile
def calculateI(nlaws,C,alpha,D,beta,Smin,Smax,Sbinlow,Sbinhigh,S0,sigma,area):
    """
    See pn_integral.pro
    This calculates I for a SINGLE bin
    #THIS DOES NOT WORK FOR NOISELESS SIMS
    """

    #I = integrate.quad(lambda x: special.jv(2.5,x), Sbinlow, Sbinhigh)
    #(I,dI) = integrate.quad(lambda S: powerLawFuncS(S,C,alpha,Smin,Smax),Sbinlow,Sbinhigh)

    # NB Answer is fairly robust to this - 1000 is a bit low so could
    # push to 5000 for production runs - think will be OK but need to
    # NB potential bias
    resolution=1000
    ds = (Smax - Smin) / resolution

    #raw_input('%f %f'%(Smin,Smax))

    I = 0.0
    # This is the integral from Smin to Smax
    # Step from Smin -> Smax
    for istep in xrange(resolution):
        Si = Smin + (istep + 0.5) * ds # + -> - because of counting from 0 not 1
        #erfs = erf((Si-Sbinlow)/(sqrt(2.0)*sigma)) - erf((Si-Sbinhigh)/(sqrt(2.0)*sigma))
        # I've added a factor of resolution here
        #I += powerLawFuncS(Si,C,alpha,Smin,Smax) * erfs * ds * resolution
        I += powerLawFuncErfsS(Si,nlaws,C,alpha,D,beta,Smin,Smax,Sbinlow,Sbinhigh,S0,\
                               sigma,area)\
          * ds * 1.0e3 # ??? XXX JZ WHAT IS THIS FACTOR OF 1000???
        #print istep,Si,I

    verbose=False
    if verbose: print I,ds,Smin,Smax,Sbinlow,Sbinhigh
    #raw_input('Press enter to continue...')

    return I

#-------------------------------------------------------------------------------

@profile
def checkConfusionNoiseOK(radioSynthOmegaSr,nlaws,C,alpha,beta,S0,gamma,S1,Smin,Smax,sigma,threshold):
    """
    Test whether confusion noise exceeds some fraction of the radio map noise
    Returns True if all is well
    """

    # Switch the prior off using checkConfusionNoise
    if not checkConfusionNoise: return True

    if confusionNoiseCheckSmax is None:
        confusionNoiseCheckSmaxUpperThresh=Smax
    else:
        confusionNoiseCheckSmaxUpperThresh=confusionNoiseCheckSmax

    confusionNoise=numpy.sqrt(calculateConfusionNoiseSq(radioSynthOmegaSr,nlaws,C,alpha,beta,S0,gamma,S1,Smin,confusionNoiseCheckSmaxUpperThresh))
    status = (confusionNoise < threshold*sigma)

    return status

#-------------------------------------------------------------------------------

@profile
def checkSminPriorOK(C,alpha,Smin,Smax,sigma,dOmega,numerical=False):

    """
    This is Ketron's equation 12
    Based on check_smin_prior.pro
    My units are in sq. deg.
    Optionally carry out the integral numerically (-> arbitrary models)
    The integral fails for Smin < 0 - but SMIN_MIN should catch this!
    """

    #checkSmin=True
    # Switch the prior off using checkSmin
    if not checkSmin: return True

    if numerical:
        JJ=integrate.quad(\
            lambda S:powerLawFuncS(S,C,alpha,Smin,Smax,SURVEY_AREA)*S**2,Smin,Smax)[0]
        #KK= integrate.quad(\
            #lambda S:powerLawFuncS(S,C,alpha,Smin,Smax,area),Smin,Smax)[0]  
        #print JJ,KK
        JJ *=  1000.0
        #print JJ
        LHS = sqrt(dOmega * JJ)
        RHS = sigma/div
        #print JJ,LHS,RHS,LHS**2,RHS**2,RHS/LHS
        #print dOmega,JJ,dOmega*JJ,(sigma/div)**2
        status = (LHS < RHS)
        if not status: print '.',
        #if not status: print C,alpha,Smin,Smax,area,LHS,RHS


        return status
    else:
        # Go the analytic route (single power law only, for now)
        pass

    #c=open('check_smin.txt','a')

    minimum_acceptable_smin=None
    try:
        minimum_acceptable_smin = ( (Smax**(3.0 + alpha)) - 
                       (((3.0 + alpha) * sigma**2) / 
                        (C * dOmega * div**2)) )**(1.0 / (3.0 + alpha))
    except ValueError:
        minimum_acceptable_smin=float('nan')

    status = (minimum_acceptable_smin <= Smin or isnan(minimum_acceptable_smin))
    if not status: print '.',

        #print '# C alpha Smin Smax sigma dOmega check_smin status'
    line= '%6.4f %6.4f %6.4f %6.4f %6.4f %e %6.4f %i\n' \
      % (C,alpha,Smin,Smax,sigma,dOmega,minimum_acceptable_smin,int(status))
    #print line

    #print minimum_acceptable_smin

      #IF finite(CHECK_SMIN) EQ 0 THEN BEGIN # infinite return 0
      #result = 1L => infinite return 1
      #RETURN, result
      #IF check_smin LE Tsmin THEN result = 1L ELSE result = 0L

      #    return check_smin,Smin

      #c.write(line)
      #c.close()

    return status
    
#-------------------------------------------------------------------------------

@profile
def myloglike(cube,ndim,nparams):
    """
    This is a Poisson lhood - eqn 6 from Ketron's paper
    compute_likelihood.pro

    There are k galaxies in the i^th bin [per sr]
    """

    #global loglike

    # This doesn't work for some reason:
    #(C,alpha,Smin,Smax) = cube
    # Try = [i for i in cube]

    C=cube[0]
    alpha=cube[1]
    Smin=cube[2]
    Smax=cube[3]
    #print C,alpha,Smin,Smax
    D=beta=S0=gamma=S1=delta=S2=-99.0
    if nlaws == 1:
        #D=-99.0
        #beta=-99.0
        #S0=-99.0
        pass
    if nlaws > 1:
        #D=cube[4]
        beta=cube[4]
        S0=cube[5]
    if nlaws > 2:
        gamma=cube[6]
        S1=cube[7]
    if nlaws > 3:
        delta=cube[8]
        S2=cube[9]

    if floatNoise:
        noise=cube[nparams-1]
    else:
        noise=sigma

    #if Smin < 100.0: print '******',Smin,Smax,C,alpha
    #print (i for i in cube) #-> MemoryError??

    #print ndim, noise
    if nlaws == 2 and not (Smin < S0 < Smax):
        print '+',
        return -1.0e99

    if nlaws == 3 and not (Smin < S0 < S1 < Smax):
        print '+',
        return -1.0e99

    if nlaws == 4 and not (Smin < S0 < S1 < S2 < Smax):
        print '+',
        return -1.0e99

    ok=checkSminPriorOK(C,alpha,Smin,Smax,noise,dOmega,numerical=numericalCheckSmin)
    if not ok:
        #print 'SminPrior not OK - cut applied'
        print '.',
        return -1.0e99

    ok2=checkConfusionNoiseOK(radioSynthOmegaSr,nlaws,C,alpha,beta,S0,gamma,S1,\
                              Smin,Smax,noise,confusionNoiseThreshold)
    if not ok2:
        print ',',
        return -1.0e99
    
    # NB ks needs to be specified if running loglike from the python prompt
    ###II=numpy.zeros(nbins-1)
    if not noisy:
        II=simBins(bins,C,alpha,Smin,Smax)
    elif noisy:
        #for ibin in xrange(nbins-1):
        #    Sbinlow=bins[ibin]
        #    Sbinhigh=bins[ibin+1]
        #    II[ibin]=calculateI(nlaws,C,alpha,D,beta,Smin,Smax,Sbinlow,Sbinhigh,S0,\
        #                        sigma,area)
        #    #print ibin,Sbinlow,Sbinhigh,II[ibin]
        #C=40.0;alpha=-1.5;Smin=1.0;Smax=20.0;noise=10.0;area=10.0
        if nlaws == 1:
            #print C,alpha,Smin,Smax
            II=calculateI3(C,alpha,Smin,Smax,area,noise=noise,\
                           dump=None,verbose=False)
        elif nlaws == 2:
            II=calculateI3(C,alpha,Smin,Smax,area,noise=noise,\
                           dump=None,verbose=False,nlaws=nlaws,\
                           D=D,beta=beta,S0=S0,gamma=gamma,S1=S1)
        elif nlaws == 3:
            II=calculateI3(C,alpha,Smin,Smax,area,noise=noise,\
                           dump=None,verbose=False,nlaws=nlaws,\
                           D=D,beta=beta,S0=S0,gamma=gamma,S1=S1)
        elif nlaws == 4:
            II=calculateI3(C,alpha,Smin,Smax,area,noise=noise,\
                           dump=None,verbose=False,nlaws=nlaws,\
                           D=D,beta=beta,S0=S0,gamma=gamma,S1=S1,\
                           delta=delta,S2=S2)


        #II *= 336273.0 / II.sum()
        #print C,alpha,Smin,Smax,area,sigma,SEED_SAMP,'xx',II.sum()

    # Two ways to calculate loglike - I've checked these are the same:
    #I=0.0
    #loglike=0.0
    #for ibin in xrange(nbins-1):
    #    I=II[ibin]
    #    k=ks[ibin]
    #    poisson = k*log(I) + k - k*log(k) - I
    #    loglike += poisson

    #if numpy.any(II<=0): print II


    #II=realiseData(nlaws,cube,area,sigma)

    # This is the same as the above and hopefully faster
    # Should iii be converted to integer??
    kk=ks[numpy.where(II > 0)]; iii=II[numpy.where(II > 0)]
    #print kk.sum(),kk.sum()*sqDeg2sr*SURVEY_AREA,kk.sum()*sqDeg2sr
    # **** This needs to be moved out of here
    #kk *= sqDeg2sr*SURVEY_AREA; iii *= sqDeg2sr*SURVEY_AREA
    #print kk,iii
    #print kk.sum(), iii.sum()
    loglike = (kk*numpy.log(iii) + kk - kk*numpy.log(kk) - iii).sum()

        #Sbinlow=bins[ibin]
        #Sbinhigh=bins[ibin+1]

        #I=calculateI(C,alpha,Smin,Smax,Sbinlow,Sbinhigh)

        #k=ksNoisy[ibin]
        #print k,I,C,alpha,Smin,Smax,Sbinlow,Sbinhigh,
        #print powerLawFuncS((Smin-Smax)/2.0,C,alpha,Smin,Smax)

        # ?? Is Smax > Smin required - ONLY I think if priors overlap??
        #if I > 0.0: #or Smin > Smax:
            #I=calculateI(C,alpha,Smin,Smax,Sbinlow,Sbinhigh)
        #else:
            #pass
        #print loglike,I,k,C,alpha,Smin,Smax,Sbinlow,Sbinhigh
        #raw_input('Press enter to continue...')

    return loglike


#-------------------------------------------------------------------------------

def readCat(cat,form=0):

    """
    Wrapper for catalogue read (to allow for varying formats)
    Array returned is ngals x [fluxes,noises]
    """
    global perGalData,perGalNgals
    
    if form==0:
        datax=numpy.genfromtxt(cat)
        ngals=numpy.shape(datax)[0]
        data=numpy.nan*numpy.ones(ngals,2)
        data[:,0]=datax[:,1]

    elif form==1:
        #cat='/Users/jtlz2/video/stacking/production4_q/pixels_Mz/pixels_all_z_0.5_1.0_Ms_9.0_9.5_Kabs_-100.0_0.0_noiseclip_1234.dat'
        datax=numpy.genfromtxt(cat)
        ngals=numpy.shape(datax)[0]
        noiseCol=10; fluxCol=12
        data=numpy.nan*numpy.ones((ngals,2))
        data[:,0]=datax[:,fluxCol]
        data[:,1]=datax[:,noiseCol]

    perGalData=data; perGalNgals=ngals
        
    return perGalData

#-------------------------------------------------------------------------------

@profile
def myprior2(cube,ndim,nparams):
    """
    """
    cube[0]=pri.UniformPrior(cube[0],C_MIN,C_MAX) # C / /Jy /deg^2
    cube[1]=pri.UniformPrior(cube[1],ALPHA_MIN,ALPHA_MAX) # alpha
    cube[2]=pri.UniformPrior(cube[2],SMIN_MIN,SMIN_MAX) # Smin / uJy
    cube[3]=pri.UniformPrior(cube[3],SMAX_MIN,SMAX_MAX) # Smax / uJy

#-------------------------------------------------------------------------------

@profile
def myloglike2(cube,ndim,nparams):

    """
    10.2.14
    First go at per-galaxy rms loglike function
    ***** DOES THIS NEED AN ERROR FUNCTION??
    """

    C=cube[0]
    alpha=cube[1]
    Smin=cube[2]
    Smax=cube[3]

    loglike=0.0
    # Loop over jgals
    for jgal in range(perGalNgals):
        jflux,jnoise = perGalData[jgal,:]
        #print Integrand2(,jflux,jnoise,C,alpha,Smin,Smax,1.0)
        Lj = integrate.quad(\
            lambda S:Integrand2(S,jflux,jnoise,C,alpha,Smin,Smax,1.0),Smin,Smax)[0]
            #print jgal,jflux,jnoise,Lj
        if Lj > 0.0:
            loglike += log(Lj)

    #loglike1 = loglike
    # ? Replace explicit loop with numpy loop (faster..?)
    #funcc=lambda perGalData,C,alpha,Smin,Smax: \
    #  IntegrateIntegrand2(perGalData[:,0],perGalData[:,1],C,alpha,Smin,Smax)
    #loglike=numpy.log(funcc(perGalData,C,alpha,Smin,Smax)).sum()
    #loglike2=loglike
    #print loglike1,loglike2

    # This is the exponential prefactor term, A (equation 11) in draft
    A = integrate.quad(\
        lambda S:powerLawFuncS(S,C,alpha,Smin,Smax,1.0),Smin,Smax)[0]
    loglike -= A

    return loglike

#-------------------------------------------------------------------------------

@profile
def IntegrateIntegrand2(jflux,jnoise,C,alpha,Smin,Smax):
    """
    """

    Lj = integrate.quad(\
        lambda S:Integrand2(S,jflux,jnoise,C,alpha,Smin,Smax,1.0),Smin,Smax)[0]

    return Lj

#-------------------------------------------------------------------------------

@profile
def Integrand2(S,jflux,jnoise,C,alpha,Smin,Smax,area):
    """
    Calculate I2 integrand for the per-galaxy rms case, for the vector
        of all galaxies
    I2 is the value of equation 10 in the draft manuscript
    """
    arg = -(S-jflux)**2 / (2.0*jnoise*jnoise)
    prefactor = 1.0/(sqrt(2.0*pi*jnoise*jnoise))
    Gflux = prefactor*exp(arg)

    integrand=Gflux*powerLawFuncS(S,C,alpha,Smin,Smax,1.0)

    return integrand

    
#-------------------------------------------------------------------------------

@profile
def func(*args): 
   print args

@profile
def dump_wrapper2(nsamples, nlive, n, postdist):
   print 'python dumper callback called!'
   sys.stdout.flush()
   print nsamples, nlive, n
   sys.stdout.flush()

@profile
def dump_wrapper(nsamples, nlive, n, postdist):
    #print nsamples, nlive, n, postdist
   arr_type = ((ctypes.c_float * nsamples) * n)
   values = arr_type(postdist)
   print nsamples, nlive, n, values
   #return

#-------------------------------------------------------------------------------



@profile
def initPlot(n_params):
    global plothandle
    plt.ion()
    plt.figure()
    if tellthetruth:
        thetruth = \
          {parameters[0]:40.0,parameters[1]:-1.5,parameters[2]:1.0,parameters[3]:20.0}
    else:
        thetruth=None
          
    if model == 1 or model == 2:
        plt.xlabel(parameters[0])
        plt.ylabel(parameters[1])
        plt.xlim(1.0,100.0)
        plt.ylim(-2.5,-0.1)

    elif model == 3:
        plt.xlabel(parameters[2])
        plt.ylabel(parameters[3]) 
        plt.xlim(0.0,5.0)
        plt.ylim(1.0*sigma,5.0*sigma)
    plothandle, = plt.plot([], [],'.',markersize=2)
    plt.draw()

    if (model == 1 or model ==2) and thetruth is not None:
        plt.plot(thetruth[parameters[0]],thetruth[parameters[1]],'+',color='red')
        plt.draw()
        
    return plothandle

@profile
def updatePlot(plothandle,cube):
    """
    http://stackoverflow.com/questions/10944621/dynamically-updating-plot-in-matplotlib
    """
    if model == 1 or model == 2:
        plothandle.set_xdata(numpy.append(plothandle.get_xdata(), cube[0]))
        plothandle.set_ydata(numpy.append(plothandle.get_ydata(), cube[1]))
    elif model == 3:
        plothandle.set_xdata(numpy.append(plothandle.get_xdata(), cube[2]))
        plothandle.set_ydata(numpy.append(plothandle.get_ydata(), cube[3]))
    plt.draw()
        #t=threading.Timer(2,plt.draw())
    return

#def updatePlot(plothandle,cube):
#    """
#    http://stackoverflow.com/questions/10944621/dynamically-updating-plot-in-matplotlib
#    """
#    #print cube[0],cube[1],cube[2],cube[3]
#    if model == 1 or model == 2:
#        plothandle.set_xdata(numpy.append(plothandle.get_xdata(), cube[0]))
#        plothandle.set_ydata(numpy.append(plothandle.get_ydata(), cube[1]))
#    elif model == 3:
#        plothandle.set_xdata(numpy.append(plothandle.get_xdata(), cube[2]))
#        plothandle.set_ydata(numpy.append(plothandle.get_ydata(), cube[3]))
#    plt.draw()
#    #t=threading.Timer(2,plt.draw())
#    return

@profile
def savePlot(outf='test.png'):
    """
    
    """
    plt.savefig(outf)
    #plt.close()
    return

#-------------------------------------------------------------------------------


@profile
def main():
    """
    model 0: vary no parameters
    model 1: vary all 4 parameters
    model 2: vary C and alpha only
    model 3: vary Smin and Smax only
    """
    global ks,noisy,sigma,area,model
    global parameters
    global verbose
    global loud,tellthetruth
    global nlaws
    global triangle
    global numericalCheckSmin

    global perGalData,perGalNgals
    
    #verbose=False

    #if len(sys.argv) < 1:
    #    print 'usage:'
    #    #print '       ./lumfunc.py NOISY  SIGMA  AREA  MODEL  LOUD  TRUTH  SEED NLAWS'
    #    #print '                    [True|False] [uJy] [sq. deg.]  [0|1|2] [True|False]'
    #    #print '                    [True|False] (suppressible by LOUD) [1|2]'
    #    #print '                    [-1 for clock]'
    #    print
    #    print 'with MPI:'
    #    print '       mpirun -n 2 ./lumfunc.py settings.py'
    #    print
    #    print 'Models:'
    #    print '       0  vary no parameters'
    #    print '       1  vary all 4 parameters'
    #    print '       2  vary C and alpha only'
    #    print '       3  vary Smin and Smax only'
    #    print
    #    sys.exit(0)

    #param_file=sys.argv[-1]
        
    #noisy=parseBoolean(sys.argv[1])
    #sigma=float(sys.argv[2])
    #area=float(sys.argv[3])
    #model=int(sys.argv[4])
    #loud=parseBoolean(sys.argv[5])
    #tellthetruth=parseBoolean(sys.argv[6])
    #seed=int(sys.argv[7])
    #nlaws=int(sys.argv[8])

    # Set up MPI
    world=MPI.COMM_WORLD
    rank=world.rank
    size=world.size
    master = rank==0
    #master=True

    if master:
        set_module=importlib.import_module(setf)
        globals().update(set_module.__dict__)

    note='MPI processors checked in: rank/size = (%i/%i)' % (rank,size)

    #outdir='chains'
    #if not os.path.exists(outdir): os.mkdir(outdir)
    # This is because each MPI process will try to create outdir - how
    # to avoid since master-only boolean didn't seem to work?
    if master:
        try:
            os.mkdir(outdir)
        except OSError:
            pass

        #logfile='README'
        logf=os.path.join(outdir,logfile)
        if master and os.path.exists(logf): os.remove(logf)
        log=open(logf,'w')
        remark(log,note)

    # Wait here after check-in...
    world.Barrier()
    if master: print 'All %i processors checked in...' % size

    # Broadcast global settings variables
    if master:
        set_dict = set_module.__dict__
    else:
        set_dict = None
        #for proc in xrange(size):
        #    print proc
        #    world.send(set_module.__dict__, dest=proc, tag=11)
        #print 'All %i processors transmitted OK...' % size
    #print set_dict
    #else:
        #set_dict = world.recv(source=0, tag=11)
    #if master:
    set_dict = world.bcast(set_dict,root=0)

    if not master:
        globals().update(set_dict)
        #print globals()

    # Wait here after broadcast...
    world.Barrier()
    if master: print 'All %i processors received OK...' % size

    #if master:
    #    comment=''
    #    comment=raw_input('Run comment: ')
    #    print comment

    # This is because everything is referenced to 10 sq. deg. at the moment
    #area *= 1.0/REF_AREA
###    multArea = area/REF_AREA

    # Set some MultiNEST parameters here
    #n_live_points=1000
    #multimodal=False
    #max_modes=1

    # Switch to INS
    #do_INS=True
    #n_live_points=500

    #noisy=True
    ks=numpy.zeros(numpy.shape(ksRaw))
    #noisy=True
    if noisy:
        # Switched off multArea
        #ks=ksNoisy * multArea / multArea
        ks=ksNoisy
        #sigma=1.0e-6 # Jy
        #sigma=10.0 # uJy
    else:
        #ks=ksRaw * multArea
        ks=ksRaw
        sigma=1.0e-90

    # Write settings variables to file
    if master:
        variablesf=os.path.join(outdir,variablesfile)
        #if os.path.exists(variablesf): os.remove(variablesf)
        dump_variable_values(set_module,variablesf,verbose=False)

    if master:
        print
        startTime = time.strftime('%X %x %Z')
        note='Time now is %s' % startTime
        remark(log,note)

        shutil.copy(param_file,outdir)
        note='Settings file: %s' % param_file
        remark(log,note)
        shutil.copy(datafile,outdir)
        note='Data file: %s' % datafile
        remark(log,note)

        # This is to allow import of settings from outdir
        # i.e. from outdir import * [or whatever]
        init_file='__init__.py'
        initf=os.path.join(outdir,init_file)
        touch(initf)

        note='Bins taken from %s' % datafile
        remark(log,note)
        note='# Bin occupancies [i uJy uJy field^-1]:'
        remark(log,note)
        for ibin in xrange(nbins-1):
            try:
                line='%i %f %f %f'%(ibin+1,bins[ibin],bins[ibin+1],ks[ibin])
            except IndexError:
                print "Probably your binstyle doesn't match the datafile bins"
                sys.exit(0)
            remark(log,line)


    #if nlaws == 1:
    #    n_params = 4
    #elif nlaws == 2:
    #    n_params = 6
    #elif nlaws == 3:
    #    n_params = 8
    n_params = 2*(nlaws+1)

    if floatNoise: n_params +=1
    
    # run MultiNest
    if master:
        t0 = time.time()

    # mode_tolerance=-1e90 is required as a bugfix (in earlier
    # versions of PyMultiNest)
        if loud:
            #        plothandle=initPlot(n_params)
            plothandle=initPlot(n_params)
            #print 'init', plothandle, plothandle.get_xdata(),plothandle.get_ydata()

    if perGalMode:
        perGalData=readCat(perGalCat,form=perGalCatForm)
        try:
            pymultinest.run(myloglike2,myprior2,n_params,resume=RESUME,verbose=True,\
                            multimodal=multimodal,max_modes=max_modes,write_output=True,\
                            n_live_points=n_live_points,\
                            evidence_tolerance=evidence_tolerance,\
                            mode_tolerance=-1e90,seed=SEED_SAMP,max_iter=max_iter,\
                            importance_nested_sampling=do_INS,\
                            outputfiles_basename=os.path.join(outdir,outstem),\
                            init_MPI=False)
        except:
            return 1
        
    else:
        try:
            #dump_callback=updatePlot(hl,(new_x,new_y))
            #progress = pymultinest.ProgressPlotter(n_params = n_params); progress.start()
            # NB MPI is already init'ed by mpi4py (crashes otherwise)
            pymultinest.run(myloglike,myprior,n_params,resume=RESUME,verbose=True,\
                            multimodal=multimodal,max_modes=max_modes,write_output=True,\
                            n_live_points=n_live_points,\
                            evidence_tolerance=evidence_tolerance,\
                            mode_tolerance=-1e90,seed=SEED_SAMP,max_iter=max_iter,\
                            importance_nested_sampling=do_INS,\
                            outputfiles_basename=os.path.join(outdir,outstem),\
                            init_MPI=False)#,sampling_efficiency='parameter')#,
            #        dump_callback=dump_wrapper2,n_iter_before_update=2)
            #new_x,new_y=-1.5,20.0
            #threading.Timer(2, updatePlot(hl,(new_x,new_y))).start()
            #progress.stop()
        except:
            return 1

    if master:
        if loud and (model == 1 or model == 2):
            outf=os.path.join(outdir,'%s-%s.png' % (parameters[0],parameters[1]))
            savePlot(outf)

        print '# Bin occupancies:'
        for ibin in xrange(nbins-1):
            print ibin+1,bins[ibin],bins[ibin+1],ks[ibin]

        t1 = time.time()
        dt=t1-t0

        note='Time then was %s' % startTime
        remark(log,note)
        stopTime=time.strftime('%X %x %Z')
        note='Time now is %s' % stopTime
        remark(log,note)
        note='Execution took %6.4f sec (~ %i min) with %i cores' % \
          (dt,int(round(dt/60.0)),size)
        remark(log,note)
        note='Arguments: %s' % ' '.join(sys.argv)
        remark(log,note)
        
        note='INS   = %s' % do_INS
        remark(log,note)
        note='nlive = %i' % n_live_points
        remark(log,note)
        note='Run comment: %s' % comment
        remark(log,note)

        note='Now execute:'
        remark(log,note)
        note='import pylab; from utils import *; import contour_plot'
        remark(log,note)
        note='from %s import settings' % outdir
        remark(log,note)
        #note="contour_plot.contourTri(pylab.loadtxt('%(od)s/%(os)spost_equal_weights.dat'),line=True,outfile='%(od)s/%(tri)s',col=('red','blue'),labels=settings.parameters,ranges=settings.plotRanges,truth=settings.plotTruth,reconstruct=(medianArray(settings.bins),lumfunc.ksRaw),autoscale=False,title='%(od)s')" \
        #% {'od':outdir,'os':outstem,'tri':triangle}
        note="contour_plot.contourTri(pylab.loadtxt('%(od)s/%(os)spost_equal_weights.dat'),line=True,outfile='%(od)s/%(tri)s',col=('red','blue'),labels=settings.parameters,ranges=settings.plotRanges,truth=settings.plotTruth,autoscale=False,title='%(od)s')" \
        % {'od':outdir,'os':outstem,'tri':triangle}

        remark(log,note)
        note='or\n./plot.py %s' % outdir
        remark(log,note)
        note='and\n./reconstruct.py %s' % outdir
        remark(log,note)

        log.close()

        # Copy the stats file so it's legible on my iPhone, Google, email etc.
        stats_dotdat= '%(od)s/%(os)sstats.dat' % {'od':outdir,'os':outstem}
        stats_dottxt= '%(od)s/%(os)sstats.txt' % {'od':outdir,'os':outstem}
        shutil.copy(stats_dotdat,stats_dottxt)
        
        # Finally, attempt to produce the triangle plot
        #pew=pylab.loadtxt('%s/%spost_equal_weights.dat'%(outdir,outstem))
        #triangle=os.path.join(outdir,triangle)
        #contour_plot.contourTri(pew,line=True,outfile=triangle,col=('red','blue'),\
        #                        labels=parameters,ranges=plotRanges,truth=plotTruth,\
        #                        autoscale=False,title=outdir)

        
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)


