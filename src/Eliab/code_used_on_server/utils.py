"""
Collection of my utilities (pertaining to lumfunc project)
Also includes some useful constants
"""

import os,sys
import shelve
import numpy
import scipy
from scipy import stats
from scipy.interpolate import interp1d
from scipy import ndimage
from profile_support import profile
import pymultinest

#-------------------------------------------------------------------------------

# Constants
from math import sqrt,pi,log
from scipy.constants import c as clight
from scipy.constants import k as kb
muJy2Jy=1.0e-6
sqDeg2sr=4.0*pi*pi/129600.0
sqrtTwo=sqrt(2.0)
Jy2muJy=1.0e6
beamFac=pi/(4.0*log(2.0))

#-------------------------------------------------------------------------------

def calculateBeamsPerSource(surveyArea_sqDeg,radioSynthBeamFWHM_arcsec,numSources):
    """
    surveyArea_sqDeg             Survey area in sq. deg.
    radioSynthBeamFWHM_arcsec    FWHM of synthesized beam in arcsec
    numSources                   Number of sources in surveyArea_sqDeg
    
    e.g. calculateBeamsPerSource(1.979178e-02,2.5,2974)
    """
    surveyArea_sr=surveyArea_sqDeg*sqDeg2sr
    radioSynthOmega_sr=sqDeg2sr*beamFac*(radioSynthBeamFWHM_arcsec/3600.0)**2
    numBeamAreas=surveyArea_sr/radioSynthOmega_sr
    numBeamAreasPerSource=numBeamAreas/float(numSources)
    sqrtNumBeamAreasPerSource=sqrt(numBeamAreasPerSource)
    
    print '#-------------------------------------------'
    print 'Survey area      / sq. deg. %4.2f' % surveyArea_sqDeg
    print '                 / sr       %e' % surveyArea_sr
    print 'Synth. beam FWHM / arcsec   %4.2f' % radioSynthBeamFWHM_arcsec
    print '            area / sr       %e' % radioSynthOmega_sr
    print 'Number of beam areas        %6.2f' % numBeamAreas
    print 'Number of sources           %i' % numSources
    print 'Number of beam areas/source %4.2f' % numBeamAreasPerSource
    print 'Number of sources/beam area %4.2f' % (1.0/numBeamAreasPerSource)
    print '       i.e. -->  %4.2f x %4.2f beams/source' % \
      (sqrtNumBeamAreasPerSource,sqrtNumBeamAreasPerSource)
    print '                (%4.2f x %4.2f sources/beam)' % \
      (1.0/sqrtNumBeamAreasPerSource,1.0/sqrtNumBeamAreasPerSource)
    print '#-------------------------------------------'

    return

#-------------------------------------------------------------------------------

@profile
def touch(fname):
    cmd='touch %s' %fname
    os.system(cmd)
    return

#-------------------------------------------------------------------------------

@profile
def touch2(fname, times=None):
    with file(fname, 'a'):
        os.utime(fname, times)
    return

#-------------------------------------------------------------------------------

def block_mean(ar, fact):
    """
    See http://stackoverflow.com/questions/18666014/downsample-array-in-python

    Run as e.g.
    ar = np.random.rand(20000).reshape((100, 200))
    block_mean(ar, 5).shape  # (20, 40)
    """
    assert isinstance(fact, int), type(fact)
    sx, sy = ar.shape
    X, Y = numpy.ogrid[0:sx, 0:sy]
    regions = sy/fact * (X/fact) + Y/fact
    res = ndimage.mean(ar, labels=regions, index=numpy.arange(regions.max() + 1))
    res.shape = (sx/fact, sy/fact)
    return res

#-------------------------------------------------------------------------------

@profile
def strictly_increasing(L):
    """http://stackoverflow.com/questions/4983258/python-how-to-check-list-monotonicity
    """
    return all(x<y for x, y in zip(L, L[1:]))

#-------------------------------------------------------------------------------

@profile
def poissonLhoodMulti2(dataObject,realisationObject,silent=True):
    #print dataObject.zsliceData
    #print realisationObject
    loglike=0.0
    #sys.exit()
    # Loop over z slices
    for r,z in enumerate(sorted(dataObject.zsliceData.keys())):
        #zbins=dataObject.zsliceData[z][1]
        data=dataObject.zsliceData[z][0]
        realisation=realisationObject[z]
        if not silent:
            for i in range(len(data)):
                print i,data[i],realisation[i]
        kk=data[numpy.where(realisation > 0)];
        iii=realisation[numpy.where(realisation > 0)]
        #loglike += (kk*numpy.log(iii) + kk - kk*numpy.log(kk) - iii).sum()
	loglike += (kk*numpy.log(iii) + kk - (kk + 0.5)*numpy.log(kk) - 0.5*numpy.log(2*numpy.pi)  - iii).sum()

    return loglike

#-------------------------------------------------------------------------------

@profile
def poissonLhoodMulti(data,realisation,silent=True,fractions=None,redshifts=None):
    loglike=0.0
    if fractions is not None:
        for j in range(fractions.size):
            loglike += poissonLhood(data[:,j],realisation[:,j],silent=silent)
        loglike += data[:,0].size * numpy.log(fractions).sum()
    elif redshifts is not None:
        print redshifts
    return loglike

#-------------------------------------------------------------------------------

@profile
def poissonLhood(data,realisation,silent=False):
    if not silent:
        for i in range(len(data)):
            print i,data[i],realisation[i]
    kk=data[numpy.where(realisation > 0)];
    iii=realisation[numpy.where(realisation > 0)]
    loglike = (kk*numpy.log(iii) + kk - (kk + 0.5)*numpy.log(kk) - 0.5*numpy.log(2*numpy.pi)  - iii).sum()
    #loglike = (kk*numpy.log(iii) + kk - kk*numpy.log(kk) - iii).sum()
    print iii,kk
    return loglike

#-------------------------------------------------------------------------------

def interpola(v, x, y,kind='linear'):
    if v <= x[0]:
        return y[0]+(y[1]-y[0])/(x[1]-x[0])*(v-x[0])
    elif v >= x[-1]:
        return y[-2]+(y[-1]-y[-2])/(x[-1]-x[-2])*(v-x[-2])
    else:
        f = interp1d(x, y, kind=kind) 
        return f(v)

#-------------------------------------------------------------------------------

def extrap1d(interpolator):
    """
    http://stackoverflow.com/questions/2745329/how-to-make-scipy-interpolate-give-an-extrapolated-result-beyond-the-input-range
    """
    xs = interpolator.x
    ys = interpolator.y

    def pointwise(x):
        if x < xs[0]:
            return ys[0]+(x-xs[0])*(ys[1]-ys[0])/(xs[1]-xs[0])
        elif x > xs[-1]:
            return ys[-1]+(x-xs[-1])*(ys[-1]-ys[-2])/(xs[-1]-xs[-2])
        else:
            return interpolator(x)

    def ufunclike(xs):
        return scipy.array(map(pointwise, scipy.array(xs)))

    return ufunclike

#-------------------------------------------------------------------------------


def interpol(func,x):
    """
    Following Russell's prescription [email of 24.1.14]
    Return an interpolation object useable for a mock generation for
        an arbitrary function/lambda

    Run as e.g.
    Ss=numpy.linspace(1.0,100.0,100)
    y=numpy.array([S**2 for S in Ss])
    lookup=buildCDF(lambda S:S**2,Ss)
    Ss_fine=numpy.linspace(1.0,100.0,1000)
    y_fine=lookup(Ss_fine)
    """

    y = numpy.array([func(ix) for ix in x])
    f = interp1d(x,y)
    return f

#-------------------------------------------------------------------------------

def buildCDF(values):
    """
    Given an array, accumulate it and normalize to the interval U[0,1]
    """
    return (values.cumsum()-values.cumsum()[0])/values.cumsum()[-1]

#-------------------------------------------------------------------------------


def matchit(x,y,shape):
    """
    http://stackoverflow.com/questions/8251541/numpy-for-every-element-in-one-array-find-the-index-in-another-array
    """
    import numpy as np
    # Every element of y must be in x
    assert np.all(np.intersect1d(x,y) == np.sort(y))
    xsorted = np.argsort(x)
    ypos = np.searchsorted(x[xsorted], y)
    indices = xsorted[ypos]
    vals=np.unravel_index(indices,shape)
    coords= [(row, ix) for row, ix in zip(vals[0],vals[1])]
    return coords

def test_matchit():
    import numpy
    numpy.random.seed(seed=1234)
    nside=4
    u=numpy.random.uniform(0,10,nside*nside).reshape(nside,nside)
    print u
    xy=numpy.random.choice(u.flatten(),size=5,replace=False)
    print xy
    coords=matchit(u.flatten(),xy,u.shape)
    print coords
    tellback=[u[i,j] for i,j in coords]
    assert(numpy.allclose(tellback,xy))
    print xy
    print tellback
    return

#-------------------------------------------------------------------------------

def p(x,xtrue):
    """
    Calculate percentage difference between two values
    """
    return 100.0*(x-xtrue)/xtrue

#-------------------------------------------------------------------------------

def gaussian(x, mu, sig,norm=True):
    gauss=numpy.exp(-0.5*((x-mu)/sig)**2)
    if norm:
        gauss *= 1.0/(sqrtTwo*sqrt(pi)*sig)
    return gauss

#-------------------------------------------------------------------------------

def sum2DGaussian(dataSlice,sigma,amplitude,x0,y0):

    """
    """
    from math import pi,sqrt,exp,log

    xlen=len(dataSlice[:,0])
    ylen=len(dataSlice[0,:])

    #print dataSlice
    gaussianSum=0.0
    for x in range(xlen):
        for y in range(ylen):
            r=sqrt((x-0.5*xlen)**2+(y-0.5*ylen)**2)
            if r > x/2.0 or r > y/2.0: continue
            amplitude = 1.0#/(pi*sigma*sigma)
            arg = -(0.5*r**2/sigma**2)
            gaussian = dataSlice[x,y] * amplitude * exp(arg)
            #print x,y,r,gaussian,dataSlice[x,y]
            gaussianSum+=gaussian

    return gaussianSum / (pi/(4.0*log(2.0))) / 4.0 / 2.0

#-------------------------------------------------------------------------------

@profile
def bilinear_interpolation(x, y, points):
    '''
    http://stackoverflow.com/questions/8661537/how-to-perform-bilinear-interpolation-in-python
    
    Interpolate (x,y) from values associated with four points.

    The four points are a list of four triplets:  (x, y, value).
    The four points can be in any order.  They should form a rectangle.

        >>> bilinear_interpolation(12, 5.5,
        ...                        [(10, 4, 100),
        ...                         (20, 4, 200),
        ...                         (10, 6, 150),
        ...                         (20, 6, 300)])
        165.0

    '''
    # See formula at:  http://en.wikipedia.org/wiki/Bilinear_interpolation

    points = sorted(points)               # order points by x, then by y
    (x1, y1, q11), (_x1, y2, q12), (x2, _y1, q21), (_x2, _y2, q22) = points

    if x1 != _x1 or x2 != _x2 or y1 != _y1 or y2 != _y2:
        raise ValueError('points do not form a rectangle')
    if not x1 <= x <= x2 or not y1 <= y <= y2:
        raise ValueError('(x, y) not within the rectangle')

    return (q11 * (x2 - x) * (y2 - y) +
            q21 * (x - x1) * (y2 - y) +
            q12 * (x2 - x) * (y - y1) +
            q22 * (x - x1) * (y - y1)
           ) / ((x2 - x1) * (y2 - y1) + 0.0)

#-------------------------------------------------------------------------------

        
@profile
def dump_variable_values(module,moduleFile,verbose=False):
    """
    http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order
    """

    # Retain only last (active) occurrence of any repeated variable
    mylist=[x for i,x in enumerate(dir(module)) if x not in dir(module)[i+1:]]

    f=open(moduleFile,'w')
    
    for i in range(len(dir(module))):
        line=['%s %s'%(x,getattr(module,x)) for x in dir(module)][i]
        if verbose: print line
        f.write('%s\n'%line)

    f.close()

    return

#-------------------------------------------------------------------------------

@profile
def peak_confidence(vector,bins=None):
    """
    For a vector, bin data and find peak position
    See http://stackoverflow.com/questions/18818905/find-the-x-value-corresponding-to-a-histogram-max
    """

    if bins is None: bins=100

    n, b = numpy.histogram(vector,bins=bins)
    bin_max = numpy.where(n == n.max())
    bin_max_posn_lower=bin_max[0][0]
    peak=(b[bin_max_posn_lower]+b[bin_max_posn_lower+1])/2.0
    
    return peak

#-------------------------------------------------------------------------------

@profile
def calculate_confidence2(vector,value_central=None,alpha=0.68,ret_all=True,\
                          truncate_edges=False):
    """
    For a given central value (could be median),
    return error bars (optionally) corrected to avoid touching the edges
    value_central is the median unless otherwise supplied (e.g. ymap[ibin])
    """
    k='strict'
    if value_central is None:
        percentile_central=50.0 # median
        value_central=stats.scoreatpercentile(vector,percentile_central)
    else:
        percentile_central=stats.percentileofscore(vector,value_central,kind=k)

    pcl=percentile_low=percentile_central-(100.0*alpha/2.0)
    pch=percentile_high=percentile_central+(100.0*alpha/2.0)

    # Correct the confidence region to avoid touching the edges
    if truncate_edges:
        if percentile_high > 100.0:
            dpc_high=percentile_high-100.0
            percentile_low -= dpc_high
            percentile_high=100.0
        elif percentile_low < 0.0:
            dpc_low=0.0-percentile_low
            percentile_low=0.0
            percentile_high += dpc_low
    #print 'www',pcl,percentile_low,pch,percentile_high,percentile_central,value_central

    #assert ((percentile_high <= 100.0) and (percentile_low >= 0.0)), '***cannot compute.... %f %f'%(percentile_low,percentile_high)

    err_low  = value_central - stats.scoreatpercentile(vector,percentile_low)
    err_high = stats.scoreatpercentile(vector,percentile_high) - value_central
    
    if ret_all:
        return value_central,err_low,err_high,\
          stats.scoreatpercentile(vector,percentile_low),\
          stats.scoreatpercentile(vector,percentile_high)
    else:
        return value_central,err_low,err_high

#-------------------------------------------------------------------------------

@profile
def calculate_confidence(vector,alpha=0.68,ret_all=False):
    """
    from stacker.py (modified)
    """

    percentile_median=50.0
    percentile_low=percentile_median-(100.0*alpha/2.0)
    percentile_high=percentile_median+(100.0*alpha/2.0)

    median   = stats.scoreatpercentile(vector,percentile_median)
    err_low  = median - stats.scoreatpercentile(vector,percentile_low)
    err_high = stats.scoreatpercentile(vector,percentile_high) - median

    if ret_all:
        return median,err_low,err_high,\
          stats.scoreatpercentile(vector,percentile_low),\
          stats.scoreatpercentile(vector,percentile_high)
    else:
        return median,err_low,err_high

#-------------------------------------------------------------------------------

@profile
def mean_confidence_interval(data, confidence=0.95):
    """
    from http://stackoverflow.com/questions/15033511/compute-a-confidence-interval-from-sample-data
    """
    a = 1.0*numpy.array(data)
    n = len(a)
    m, se = numpy.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

#-------------------------------------------------------------------------------

@profile
def randomp(power,n,range_x=[5,100],seed=None,dump=None):

    """
    This is my version of
    http://idlastro.gsfc.nasa.gov/ftp/pro/math/randomp.pro
    since scipy.stats.powerlaw doesn't seem to do the job
    Use 'dump' to specify a dumpfile
    """

    import numpy

    if len(range_x) < 1: range_x=[5,100]
    elif len(range_x) != 2:
        print '*** Error - range_x keyword must be a 2-element vector'
        return

    pow1 = power + 1.0
    lo = range_x[0] ; hi = range_x[1]
    if lo > hi:
        temp=lo
        lo=hi
        hi=temp

    if seed is not None:
        numpy.random.seed(seed=seed)
    r = numpy.random.rand(n)

    if power != -1.0:
        norm = 1.0/(hi**pow1 - lo**pow1)
        expo = numpy.log10(r/norm + lo**pow1)/pow1
        x = 10.0**expo
    else:
        norm = 1.0/(numpy.log(hi) - numpy.log(lo))
        x = numpy.exp(r/norm + numpy.log(lo))

    if dump is not None:
        #dumpf='fluxes.txt'
        dumpf=dump
        numpy.savetxt(dumpf,x)
        print 'Draws are in %s' % dumpf

    return x

#-------------------------------------------------------------------------------

@profile
def medianArray(bins):
    """
    I want to use map/reduce/ufuncs etc. to compute the median bins
    But for now just do this manually...
    """
    import numpy
    #bins=numpy.array([1.000,1.394,1.945,2.714,3.786,5.281,7.368,10.27,14.33,20.00])
    nbins=len(bins)
    
    bin_medians=numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        posts=(bins[ibin],bins[ibin+1])
        bin_medians[ibin]=numpy.median(posts)
        #print ibin,bins[ibin],bin_medians[ibin]

    return bin_medians

#-------------------------------------------------------------------------------

@profile
def meanArray(bins):
    """
    I want to use map/reduce/ufuncs etc. to compute the median bins
    But for now just do this manually...
    """
    import numpy
    nbins=len(bins)
    
    bin_means=numpy.zeros(nbins-1)
    for ibin in xrange(nbins-1):
        posts=(bins[ibin],bins[ibin+1])
        bin_means[ibin]=numpy.mean(posts)

    return bin_means

#-------------------------------------------------------------------------------

@profile

def runHistogramModes(globstring,index=0,nbins=1000):
    """
    utils.runHistogramModes('/Users/jtlz2/Dropbox/lumfunc/chains_140331[c-l]/1-post_equal_weights.dat',index=3,nbins=1000)
    utils.runHistogramModes('/Users/jtlz2/Dropbox/lumfunc/chains_140331c/1-post_equal_weights_cl.dat',index=1,nbins=1000)
    utils.runHistogramModes('/Users/jtlz2/Dropbox/lumfunc/chains_140328[k-t]/1-post_equal_weights_kt.dat',index=2,nbins=100)
    """
    import glob

    files=glob.glob(globstring)

    g=[]
    for f in files:
        g.append(histogramMode(f,index,nbins))

    m=numpy.mean(g); n=numpy.median(g); s=numpy.std(g)

    print 'Mean peak: %f +/- %f (median %f)' % (m,s,n)

    return
    
#-------------------------------------------------------------------------------

@profile
def histogramMode(datafile,index=0,nbins=1000,bins=None):
    """
    """

    import numpy

    data=numpy.genfromtxt(datafile)

    reversed=data[::-1,index]
    data[:,index]=reversed
    
    if index == 0: binmin=0.0; binmax=100.0
    elif index == 1: binmin=-2.5; binmax=-0.1
    elif index == 2: binmin=1.0; binmax=500.0
    elif index == 3: binmin=500.0; binmax=1000.0
    
    
    #bins=numpy.linspace(data[:,index].min(),data[:,index].max(),nbins)
    #if bins is None:
    bins=numpy.linspace(binmin,binmax,nbins)

    counts,bins=numpy.histogram(data[:,index],bins=bins)
    #print counts
    #print bins
    #print counts.max(),counts.argmax()
    peak=bins[counts.argmax()]

    return peak

#-------------------------------------------------------------------------------

@profile
def tail(thefile):
    f=open(thefile,'r')
    time.sleep(0.1)
    lines=f.readlines()
    #print lines
    f.seek(0,2)
    #line=f.read()
    f.close()
    if lines != []:
        return lines[-1]
    else:
        return None

@profile
def follow(thefile):
    """
    http://stackoverflow.com/questions/1475950/tail-f-in-python-with-no-time-sleep
    """
    thefile.seek(0,2)      # Go to the end of the file
    while True:
         line = thefile.readline()
         if not line:
             time.sleep(0.1)    # Sleep briefly
             continue
         yield line

#-------------------------------------------------------------------------------

def find_nearest(array,value):
    """
    http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """

    idx = (numpy.abs(array-value)).argmin()
    return idx,array[idx]

#-------------------------------------------------------------------------------


@profile
def parseBoolean (b):
    """
    http://codecomments.wordpress.com/2008/04/08/converting-a-string-to-a-boolean-value-in-python
    """
    # Handle case where b is already a Boolean type.
    if b == False or b == True:
        return b
    b = b.strip()
    if len(b) < 1:
        raise ValueError ('Cannot parse empty string into boolean')
    b = b[0].lower()
    if b == 't' or b == '1' or b == 'y':
        return True
    return False

#-------------------------------------------------------------------------------

@profile
def remark(fHandle,note,verbose=True):
    """
    print to screen (if verbose) and logfile
    """

    if verbose:
        print note
    fHandle.write('%s\n'%note)
    fHandle.flush()

    return

#-------------------------------------------------------------------------------

@profile
def remarks(fHandle,notes,verbose=True):
    """
    Loop over remark function
    """
    for note in notes:
        remark(fHandle,note,verbose=verbose)
    return

#-------------------------------------------------------------------------------

def reportRelativeEvidences(globList):
    """
    A more useful version of reportRelativeEvidences
    """

    print '\n-> Comparing runs %s' % (', '.join(globList))
    Hdict={}; Zdict={};dZdict={}#; shelvesDict={}
    for run in globList:
        statsf='%s/1-'%run
        nparams=numpy.genfromtxt('%spost_equal_weights.dat'%statsf).shape[-1]-1
        ana=pymultinest.analyse.Analyzer(nparams,outputfiles_basename=statsf).get_stats()
        Z=ana['global evidence']; dZ=ana['global evidence error']
        Hdict[Z]=run; Zdict[run]=Z; dZdict[run]=dZ
        #try:
        #    shelvesDict=restoreShelf('%s/shelves.txt'%run)
        #except:
        #    print 'No shelves found'

    #Identify H0
    Zmin=min(Hdict.keys())
    Zmax=max(Hdict.keys())
    H0=Hdict[Zmin]; Hmax=Hdict[Zmax]
    print
    print 'H0   = %s (Z=%f)' % (H0,Zmin)
    print 'Hmax = %s (Z=%f)' % (Hmax,Zmax)
    print

    print '# H_i Z dZ DeltaZ dDeltaZ null verdict'    
    for run in globList:
        DeltaZ=Zdict[run]-Zdict[H0]
        if False:#shelvesDict[run]['doRayleigh']:
            Ry='R'
        else:
            Ry=''
        if run==H0:
            dDeltaZ=0.0
        else:
            dDeltaZ=numpy.sqrt(dZdict[run]**2+dZdict[H0]**2)
        if run==H0:
            result='<-*- H0'
        elif run==Hmax:
            result='<-*- winner'
        else:
            result=' '
        print '%s: Z = %f +/- %6.4f | %6.4f +/- %6.4f %s %s' \
          % (run,Zdict[run],dZdict[run],DeltaZ,dDeltaZ,Ry,result)
    print

    print "% ...and here's the LaTeX:"
    print '$H_i$ & $\Delta\log_{\mathrm{e}}Z$ \\\\ % run'    
    for irun,run in enumerate(globList):
        DeltaZ=Zdict[run]-Zdict[H0]
        if run==H0:
            dDeltaZ=0.0
        else:
            dDeltaZ=numpy.sqrt(dZdict[run]**2+dZdict[H0]**2)
        print '%i & $%6.4f \pm %6.4f$ \\\\ %% %s' \
          % (irun+1,DeltaZ,dDeltaZ,run)
    print
    print '#% Happy now..?\n'

    #for run in globList:
    #    if run==H0:continue
    #    reportRelativeEvidence(H0=H0,H1=run,verbose=True)

    return
        
#-------------------------------------------------------------------------------

def reportRelativeEvidence(H0=None,H1=None,verbose=True):
    """
    Print Delta Z report for H0: Z0
                             H1: Z1
    Watch out for signs
    """

    if verbose:
        print '-> Comparing runs %s and %s' % (H0,H1)

    statsf0='%s/1-'%H0; statsf1='%s/1-'%H1
    nparams0=numpy.genfromtxt('%spost_equal_weights.dat'%statsf0).shape[-1]-1
    nparams1=numpy.genfromtxt('%spost_equal_weights.dat'%statsf1).shape[-1]-1

    ana0=pymultinest.analyse.Analyzer(nparams0,outputfiles_basename=statsf0).get_stats()
    ana1=pymultinest.analyse.Analyzer(nparams1,outputfiles_basename=statsf1).get_stats()
    Z0=ana0['global evidence']; dZ0=ana0['global evidence error']
    Z1=ana1['global evidence']; dZ1=ana1['global evidence error']

    DeltaZ10=Z1-Z0
    dDeltaZ10=numpy.sqrt(dZ1**2+dZ0**2)

    assert(Z1>Z0), 'hypotheses inverted, Z1 (%f) < Z0 (%f) - try again' % (Z1,Z0)

    if verbose:
        print
        print 'Model 0 - %s     : Z0 = %f +/- %f' %(H0,Z0,dZ0)
        print 'Model 1 - %s     : Z1 = %f +/- %f' %(H1,Z1,dZ1)
        print 'Log-Evidence in favour of Model 1 = %f +/- %f' % (DeltaZ10,dDeltaZ10)
        print '                                 at %f sigma' % abs(DeltaZ10/dDeltaZ10)
        print 'The odds ratio is %f:1' % numpy.exp(abs(DeltaZ10))
        print '        i.e. 10^{%f}:1' % numpy.log10(numpy.exp(abs(DeltaZ10)))
        print
        print 'LaTeX:'
        if numpy.log10(numpy.exp(abs(DeltaZ10))) > 1:
             print 'A & $%f \pm %f$ & $10^{%f}$:1 \\\\' \
               % (DeltaZ10,dDeltaZ10,numpy.log10(numpy.exp(abs(DeltaZ10))))
        else:
            print 'A & $%f \pm %f$ & %f:1 \\\\' \
               % (DeltaZ10,dDeltaZ10,numpy.exp(abs(DeltaZ10)))
    else:
        print H1

    return

#-------------------------------------------------------------------------------

@profile
def fetchStats(outdir,parameters,truth):
    """
    """

    print '-> analysing summary stats:'

    n_params=len(parameters)
    statsf='%s/1-'%outdir

    x=pymultinest.analyse.Analyzer(n_params,outputfiles_basename=statsf)
    y=x.get_stats()
    bf=x.get_best_fit()
    
    Zg=y['global evidence']
    stats=y['marginals']

    B=bf['log_likelihood']

    summary={}
    for ip,param in enumerate(parameters):
        s=stats[ip]
        b=bf['parameters'][ip]
        print s['1sigma']
        summary[param]=(b,s['median'],s['1sigma'][0],s['1sigma'][-1])

    # ugliest syntax ever!
    print '\n# truth param bestfit median lower upper'
    for param in parameters:
        print '%7s'%param,'%.2f'%truth[param],' '.join(['%.2f'%s for s in summary[param]])

    print '****Global log-evidence is %f' % Zg

    return summary

#-------------------------------------------------------------------------------

@profile
def printLaTeX(parameters,statsDict,dump=None):
    """
    ' \\\\\n'.join([' & '.join(map(str,line)) for line in a])
    """

    if dump is not None:
        outf='%s/params.tex'%dump
        out=open(outf,'w')

    for ip,param in enumerate(parameters):
        val=statsDict[param]
        line = '%6s & $%5.2f_{%6.3f}^{%6.3f}$ \\\\' % (param,val[0],val[2],val[3])
        if dump is not None:
            out.write('%s\n'%line)
        else:
            print line

    if dump is not None:
        out.close()
        print '\n-> writing summary stats to \input{%s}'%outf

    return

#-------------------------------------------------------------------------------

@profile
def restoreShelf(shelvef):
    """
    http://stackoverflow.com/questions/2960864/how-can-i-save-all-the-variables-in-the-current-python-session
    """
    my_shelf = shelve.open(filename)
    shelvesDict={}
    for key in my_shelf:
        shelvesDict[key]=my_shelf[key]
    my_shelf.close()
    return shelvesDict

#-------------------------------------------------------------------------------
