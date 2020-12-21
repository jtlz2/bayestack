"""
Support functions for bayestack, bayestackClasses and lumfunc

Jonathan Zwart
May 2015

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
from lumfuncUtils import dNdS_LF,LF, get_dlds, get_Vmax,get_sbins
import matplotlib.pyplot as plt
#from bayefit import dataModel,Lopt2Lrad
from profile_support import profile
from utils import sqDeg2sr,sqrtTwo,find_nearest,medianArray,\
                           interpol,buildCDF,Jy2muJy,interpola

if 'chains' in sys.argv[-1]:
    potential_settings=glob.glob(os.path.join(sys.argv[-1],'*settings*py'))
    #assert(len(potential_settings)==1), '***More than one potential settings file!'
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
def simulate(family,params,paramsList,bins,\
             seed=None,N=None,noise=None,output=None,\
             dump=None,version=2,verbose=False,area=None,\
             skadsf=None,pole_posns=None,simarrayf=None,\
             simdocatnoise=True):
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

    Need to add normalization capability
    
    Families:
    ========

    skads:
    -----

    r=countUtils.simulate('skads',[0.01,85.0],['S0','S1'],numpy.linspace(-60.0,100.0,26),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)
    
    ppl:
    ---

    r=countUtils.simulate('ppl',[1000.0,5.0,75.0,-1.6],['C','S0','S1','a0'],numpy.linspace(-20.0,100.0,22),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)

    r=countUtils.simulate('ppl',[1000.0,5.0,25.0,75.0,-1.6,-2.5],['C','S0','S1','S2','a0','a1'],numpy.linspace(-20.0,100.0,22),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)

    r=countUtils.simulate('ppl',[1000.0,5.0,25.0,40.0,75.0,-1.6,-2.5,-1.0],['C','S0','S1','S2','S3','a0','a1','a2'],numpy.linspace(-20.0,100.0,22),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)

    r=countUtils.simulate('ppl',[1000.0,5.0,25.0,40.0,75.0,90.0,-1.6,-2.5,-1.0,2.0],['C','S0','S1','S2','S3','S4','a0','a1','a2','a3'],numpy.linspace(-20.0,100.0,22),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)

    poly:
    ----

    r=countUtils.simulate('poly',[5.0,75.0,1.0],['S0','S1','p0'],numpy.linspace(-20.0,100.0,22),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)

    r=countUtils.simulate('poly',[5.0,75.0,1.0,-1.0],['S0','S1','p0','p1'],numpy.linspace(-20.0,100.0,22),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)
    
    r=countUtils.simulate('poly',[5.0,75.0,1.0,-1.0,5.0],['S0','S1','p0','p1','p2'],numpy.linspace(-20.0,100.0,22),seed=1234,N=40000,noise=17.0,dump='R.txt',output='dummy.txt',verbose=True)
    
    bins:
    ----

    

    test:
    ----

    array:
    -----


    """

    # Initialize seed for variates AND any noise
    if seed is not None:
        numpy.random.seed(seed=SEED_SIM)

    if family=='ppl':
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

        function = lambda S:powerLawFuncWrap(nlaws,S,C,alpha,-99.0,beta,\
                                      Smin/1e6,Smax/1e6,S0/1e6,gamma,S1/1e6,delta,S2/1e6,1.0)

    elif family=='test':
        Smin=params[paramsList.index('S0')]
        Smax=params[paramsList.index('S1')]
        function = lambda S:S**2

    elif family=='poly':
        Smin=params[paramsList.index('S0')]
        Smax=params[paramsList.index('S1')]
        coeffs=[params[paramsList.index(p)] for p in paramsList if p.startswith('p')]
        S_1=1.0
        function = lambda S:polyFunc(S,S_1,Smin,Smax,coeffs)

    elif family=='bins':
        Smin=params[paramsList.index('S0')]
        Smax=params[paramsList.index('S1')]
        coeffs=[params[paramsList.index(p)] for p in paramsList if p.startswith('b')]
        if pole_posns is None:
            pole_posns=numpy.logspace(numpy.log10(Smin),numpy.log10(Smax),len(coeffs)+1)
        assert(len(coeffs)==len(pole_posns)-1), '***Mismatch in number of poles!!'
        Smin=pole_posns[0]
        Smax=pole_posns[-1]
        function = lambda S:polesFunc(S,pole_posns,Smin,Smax,coeffs)

    elif family=='array':
        Smin=params[paramsList.index('S0')]
        Smax=params[paramsList.index('S1')]
        assert(simarrayf is not None), '***Need to specify an input simulation!'
        print 'Reading %s...' % simarrayf
        dataMatrix=numpy.genfromtxt(simarrayf)
        dndsInArr=dataMatrix[:,4]
        binsDogleg=numpy.concatenate((dataMatrix[:,0],[dataMatrix[-1,1]]))
        binsMedian=dataMatrix[:,2]
        assert((medianArray(binsDogleg)==binsMedian).all()), '***bin mismatch!'
        Smin=binsDogleg[0]; Smax=binsDogleg[-1]
        if not simdocatnoise:
            Smin=-5.01#-2.01 # binsMedian[0]
        print dndsInArr
        function=lambda S:arrayFunc(S,binsMedian,dndsInArr,Smin,Smax)

        #function2=lambda S:arrayFunc(S,binsMedian,dndsInArr,Smin,Smax)
        #for x in numpy.linspace(-10.0,100.0,500):
        #    print x,function(x),function2(x)
        #sys.exit(0)

    elif family=='skads':
        Smin=params[paramsList.index('S0')]
        Smax=params[paramsList.index('S1')]
        function=None
        assert(skadsf is not None), '***Need to specify input SKADS file!'
        print 'Reading %s...' % skadsf
        R=Jy2muJy*10**numpy.genfromtxt(skadsf)
        numpy.ndarray.sort(R)
        iRmin,Rmin=find_nearest(R,Smin)
        iRmax,Rmax=find_nearest(R,Smax)
        F=R[iRmin:iRmax]
        print '%i/%i sources ingested after Smin/Smax cuts' % (len(F),len(R))
        if N is not None:
            F=numpy.random.choice(F,size=N,replace=False)
        N=len(F)
        print 'NSKADS = %i' % N
        
    elif family =='Lrad':
    	Smin=params[paramsList.index('LoptMIN')]
    	Smax=params[paramsList.index('LoptMAX')]
    	A=params[paramsList.index('A')]
    	B=params[paramsList.index('B')]
    	sigma_Lrad=params[paramsList.index('sigma_Lrad')]
    	#print Loptmin,Loptmax
    	print 'Doing LF simulation'
    	inta =None
    	#intg = integrate.quad(lambda Lopt:Lopt2Lrad(Lopt,A=A,B=B,flux=False),Loptmin,Loptmax,epsabs=0.)[0]

    	function = lambda Lopt:Lopt2Lrad(Lopt,A=A,B=B,flux=False)
    elif family in ['LFsch','LFdpl']:
        redshift = 0.325
        z_min = 0.2
        z_max = 0.45
    	Lmin=params[paramsList.index('LMIN')]
    	Lmax=params[paramsList.index('LMAX')]
    	[Smin,Smax]= SMIN_SIM,SMAX_SIM
    	print Smin,Smax
    	[Smin,Smax]= get_sbins([10**Lmin,10**Lmax],redshift,dl)*1e6
    	print Smin,Smax,Lmin,Lmax
    	print 'Doing LF simulation'
    	Vmax=get_Vmax(z_min,z_max)
    	dsdl = get_dsdl(redshift,dl)
    	inta =None
    	intg = integrate.quad(lambda S:LF(S,redshift,dsdl,Vmax,dl,params=params,paramsList=paramsList,\
             inta=inta,area=area,family=family),Smin*1e-6,Smax*1e-6,epsabs=0.)[0]
        print intg*Vmax
        print Vmax
        area = N/(Vmax*intg)
        area1=area
        print N,area

    	function = lambda S:dNdS_LF(S,z_min,redshift,z_max,dl,params=params,paramsList=paramsList,\
             area=area,family=family)    	 

    if family != 'skads':
        # Set up the 'rough' array
        gridlength=10000 # Good enough to prevent bleeding at the edges
        Ss=numpy.linspace(Smin,Smax,gridlength)
        print Smin,Smax
        print 'checking for one sample'
        kl = function(20/1e6)
        print kl
        #sys.exit()
        values=numpy.array([function(ix/1e6) for ix in Ss])
        print values[:10]
        # Build the CDF
        CDF=buildCDF(values)
        plt.plot(CDF,Ss)
        #plt.xscale('log')
        #plt.yscale('log')
        plt.ylabel('Flux')
        plt.xlabel('CDF')
        plt.show()
        print CDF.max()
        # Create the interpolant object
        sampler=interp1d(CDF,Ss)
        
        plt.plot(Ss,values,'.')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Flux')
        plt.ylabel('LF')
        plt.show()
        
        x = numpy.linspace(0.,1.,10000)
        z = numpy.logspace(0,1,1000)/10.
        f = sampler(z)
        y = sampler(x)
        plt.yscale('log')
        #plt.xscale('log')
        plt.axhline(Smin)
        plt.axhline(Smax)
        plt.xlabel('R')
        plt.ylabel('Sampler(R) [flux]')
        #plt.plot(x,y)
        plt.plot(z,f)
        plt.show()
        #sys.exit()

        # Test that the sampler extrema match
        print Smin,sampler(0.0), 'you know wa mean'
        print Smax,sampler(0.99999)
#        assert(numpy.isclose(sampler(0.0),Smin)[0])
#        assert(numpy.isclose(sampler(0.99999),Smax,atol=1.0e-3)[0])

        # Draw the random deviates
        R = numpy.random.rand(N)
        print len(R)
        F=sampler(R)
        Nt = 0.
        for f in F:
         if f<1.:
          Nt+=1.
        F = F[F>1.]
        print N,Nt
        print len(F)
        Nt = len(F)
        #sys.exit()
        
        # Normalize here - this is N2C
        # EITHER N is specified explicitly
        # BOTH N2C and C2N are useful
        # Integrate the original function
        #intg = integrate.quad(lambda S:LF(S,redshift,dsdl,Vmax,dl,params=params,paramsList=paramsList,\
         #    inta=inta,area=area,family=family),Smin*1e-6,Smax*1e-6,epsabs=0.)[0]
        #print intg*Vmax
        #print Vmax
        #area = Nt/(Vmax*intg)
        #print N,area,area1
        #plt.show()
        #sys.exit()
        A=integrate.quad(function,Smin,Smax)[0]
#        print A,N
        # Bin the random samples
        bbins=numpy.linspace(Smin,Smax,100)
        E=numpy.histogram(F,bins=bbins)[0]
        # And calculate their area
        G=integrate.trapz(E,x=medianArray(bbins))
#        print G
#        print G/A
        # Gunpowder, treason and....
        if False:
            plt.xlim(0.0,100.0)
            plt.xlabel('S / $\mu$Jy')
            plt.hist(F,bins=bbins)
            plt.plot(Ss,values*G/A,'r')
            plt.savefig('N2C.pdf')
            plt.close()

    # Want: C given N, to compare to original C
    numbins=1000
    if family=='ppl':
        C_calc=N/N2C(function,F,Smin,Smax,numbins)
        #print N2C(function,F,Smin,Smax,numbins),C
        print 'For %i sources, C is %e (should be %e)' % (N,C_calc,C)
    elif family=='poly':
        C_calc=log10(N/N2C(function,F,Smin,Smax,numbins))
        print 'For %i sources, C is %e (should be %e)' % (N,C_calc,coeffs[0])

    # Dump noiseless fluxes to file
    puredumpf=dump
    idl_style=False
    numpy.savetxt(puredumpf,F)
    print 'Draws (noiseless) are in %s' % puredumpf
    writeCountsFile(output[1],bins,F,area,idl_style=idl_style,verbose=verbose)
    print output[1]
    # Now add noise if requested
    if simdocatnoise:
        numpy.random.seed(seed=SEED_SIM)
        poln=False
        if poln:
            F+=rice.rvs(F/noise,size=N)
        else:
            F+=numpy.random.normal(0.0,noise,Nt)

    # Dump noisy fluxes to file
    if dump is not None:
        noisydumpf='%s_noisy.txt' % puredumpf.split('.')[0]
        numpy.savetxt(noisydumpf,F)
        print 'Draws (noisy) are in %s' % noisydumpf
        print 'Minimum flux in catalogue = %f' % F.min()
        print 'Maximum flux in catalogue = %f' % F.max()

    # Write counts file
    print output[0]
    writeCountsFile(output[0],bins,F,area,idl_style=idl_style,verbose=verbose)
    print N,area#,area1

    return F

#-------------------------------------------------------------------------------

@profile
def N2C(function,deviates,Smin,Smax,numbins):
    """
    Since C is a function of the number of deviates N drawn from the
    function, calculate what the ratio is and return it
    """
    A=integrate.quad(function,Smin,Smax)[0]
#   print A,N
    # Bin the random samples
    bbins=numpy.linspace(Smin,Smax,numbins)
    E=numpy.histogram(deviates,bins=bbins)[0]
    # And calculate their area
    G=integrate.trapz(E,x=medianArray(bbins))

    return numbins*G/A

#-------------------------------------------------------------------------------

@profile
def writeCountsFile(output,bins,fluxes,area,idl_style=None,\
                    version=2,verbose=None,corrs=None, redshift=None):
    """
    Write an array of fluxes to a binned counts file
    """

    # Test version
    if version < 2: return '***Unsupported!!'

    # Bin up the fluxes
    counts=numpy.histogram(fluxes,bins=bins)[0]
    index = numpy.digitize(fluxes,bins = bins)
    if redshift is not None:
       red_mean = numpy.zeros(len(bins))
       for i in range(len(bins)):
          red_mean[i] = numpy.mean(redshift[index ==i])
          #print numpy.mean(redshift[index ==i])         
    N=len(fluxes)
    print '-> %i/%i objects observed in total (after binning)\n' % (counts.sum(),N)
    
    # Calculate differential counts
    idl_style=False
    print bins
    
    #dn_by_ds,dn_by_ds_eucl,dn_by_ds_errs,dn_by_ds_b,dn_by_ds_b_errs=\
    #  calculateDnByDs(1.0e-6*bins,counts,idl_style=idl_style,return_all=True)
    #dn_by_ds,dn_by_ds_eucl,dn_by_ds_errs,dn_by_ds_b,dn_by_ds_b_errs=\
    #  calculateDnByDs(1.0e-6*bins,counts,idl_style=idl_style,return_all=True,verbose=1)

    median_bins=medianArray(bins) # uJy
    NB=len(median_bins)

    # Set up the corrections
    if corrs is None:
        corrs=numpy.ones(NB)	

    if output is not None:
        outputf=output
        s=open(outputf,'w')
        if version < 2:
            header='# bin_median ksRaw ksNoisy'
        elif redshift is not None:
              header='# bin_low_uJy bin_high_uJy bin_median_uJy Ns_tot_obs Redshift'

        else:
            header='# bin_low_uJy bin_high_uJy bin_median_uJy Ns_tot_obs dnds_srm1Jym1 dnds_eucl_srm1Jy1p5 delta_dnds_eucl_lower_srm1Jy1p5 delta_dnds_eucl_upper_srm1Jy1p5 corr Ncgts_degm2 dNcgts_lower_degm2 dNcgts_upper_degm2'
        s.write('%s\n'%header)
        if verbose: print header
        for ibin in range(NB -1):
            if version < 2:
                line='%f %i %i' % (median_bins[ibin],-99.0,counts[ibin])
                
            elif redshift is not None:
                    line='%f %f %f %i %f' % (bins[ibin],bins[ibin+1],median_bins[ibin],counts[ibin],red_mean[ibin+1])
                   
           
            else:
                     
                line='%f %f %f %i %e %e %e %e %f %i %i %i' % \
                  (bins[ibin],bins[ibin+1],median_bins[ibin],round(counts[ibin]),\
                   dn_by_ds[ibin]/(sqDeg2sr*area),\
                   dn_by_ds_eucl[ibin]/(sqDeg2sr*area),\
                   dn_by_ds_errs[ibin]/(sqDeg2sr*area),\
                   dn_by_ds_errs[ibin]/(sqDeg2sr*area),\
                   corrs[ibin],\
                   round(counts[ibin:].sum()*1.00/area),\
                   round(numpy.sqrt(counts[ibin:].sum()*1.00/area)),\
                   round(numpy.sqrt(counts[ibin:].sum()*1.00/area)))
            s.write('%s\n'%line)
            if verbose: print line
        print counts.sum()
        s.close()

        print 'Look in %s' % outputf

    return

#-------------------------------------------------------------------------------

@profile
def calculateDnByDs(bins,counts,eucl=True,verbose=False,idl_style=False,
                    errors=False,bright=False,bright_errors=False,
                    return_all=False):
    """

    This function expects bins to be in Jy
    The output units are Jy^-1, or Jy^1.5 if eucl=True

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
        if verbose:
         print 'dnds',dn/ds
         print 'dnds eucl', (Smed**2.5) * dn/ds
         print 'dnds bright', (Smed**2.0) * dn/ds
         #return dn/ds, (Smed**2.5) * dn/ds,\
          #(Smed**2.5) * numpy.sqrt(dn)/ds,\
          #(Smed**2.0) * dn/ds,\
          #(Smed**2.0) * numpy.sqrt(dn)/ds
          
        return dn/ds, (Smed**2.5) * dn/ds,\
          (Smed**2.5) * numpy.sqrt(dn)/ds,\
          (Smed**2.0) * dn/ds,\
          (Smed**2.0) * numpy.sqrt(dn)/ds,\
          Smed


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

    return erfs * polyFunc(S,S_1,Smin,Smax,coeffs) * aarea

#-------------------------------------------------------------------------------

@profile
def polesFuncErfsS(S,pole_posns,coeffs,Smin,Smax,Sbinlow,Sbinhigh,ssigma,aarea):
    """
    """
    
    if S < Smin or S > Smax:
        return 0.0

    erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

    return erfs * polesFunc(S,pole_posns,Smin,Smax,coeffs) * aarea

#-------------------------------------------------------------------------------

@profile
def powerLawFuncS(S,C,alpha,Smin,Smax,area):
    """
    Ketron's equation (9)
    """

    if S <= Smin or S >= Smax:
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
def polyFunc(S,S_1,Smin,Smax,c):
    """
    """

    if S <= Smin or S >= Smax:
        return 0.0
    #assert(len(c)==3)
    exponent=0.0
    for n in range(len(c)):
        #exponent += c[n] * (numpy.log10(S/S_1)**n)
        exponent += c[n]*(S**n)
    #print c
    return exponent
    #return 10**exponent

#-------------------------------------------------------------------------------

@profile
def arrayFunc(S,binCentres,dndsInArray,Smin,Smax):
    """
    Truncate outside [Smin,Smax]
    AND
    linearly extrapolate outside bin centres
    """

    if S <= Smin or S >= Smax:
        return 0.0

    return interpola(S,binCentres,dndsInArray,kind='linear')

#-------------------------------------------------------------------------------

@profile
def polesFunc(S,pole_posns,Smin,Smax,coeffs):
    """
    """

    if S <= Smin or S >= Smax:
        return 0.0
    
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
    return coeffs[idx]#/(pole_posns[idx+1]-pole_posns[idx])

#-------------------------------------------------------------------------------

@profile
def calculateI(params,paramsList,bins=None,area=None,
                family=None,dump=None,verbose=False,model=None):

    """
    pn_integral, but for various different function families
    Flux arguments to this function are in uJy
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
            II[ibin]=integrate.quad(lambda S:polesFuncErfsS(S,pole_posns/1.0e6,coeffs,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2sr*area),pole_posns[0]/1.0e6,pole_posns[-1]/1.0e6)[0]
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
