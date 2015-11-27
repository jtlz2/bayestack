#!/usr/bin/env python

"""
This is reconstruct.py
Jonathan Zwart
May 2015

Usage:

./reconstruct.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
from scipy import stats
import pymultinest
from bayestackClasses import countModel
from utils import sqDeg2sr,fetchStats,calculate_confidence,\
     calculate_confidence2,peak_confidence,medianArray

param_file=sys.argv[-1]
settingsf='%s.bayestack_settings' % param_file

#-------------------------------------------------------------------------------

def main():

    """
    """

    # Import the settings variables
    print 'Settings file is %s' % param_file

    # Import the settings variables
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)

    # Set up the experiment
    # A first pass is needed in order to fetch the parameter mapping
    dummy=countModel(modelFamily,nlaws,settingsf,dataset,floatNoise,\
                    doPoln=doPoln,doRayleigh=doRayleigh)
    # Now we can fetch SMIN_MAP
    plotTruth=dict((name,-99.0) for name in dummy.parameters)
    SMIN_MAP=fetchStats(outdir,dummy.parameters,plotTruth)['S0'][0]

    startBin=SMIN_MAP; stopBin=5.0*SURVEY_NOISE; nbins=100 # For evaluation + plotting
    mybins=numpy.logspace(numpy.log10(startBin),numpy.log10(stopBin),nbins)
    expt=countModel(modelFamily,nlaws,settingsf,dataset,floatNoise,\
                    doPoln=doPoln,doRayleigh=doRayleigh,mybins=mybins)

    f='%spost_equal_weights.dat' % outstem
    f=os.path.join(outdir,f)

    # Load equally-weighted posterior samples
    x=numpy.genfromtxt(f)
    nsamp=x.shape[0]
    ncols=x.shape[1] # The fifth [seventh] column is the posterior value
    # There must be a better way, but:
    z=numpy.zeros((nsamp,ncols-1+expt.nbins))
    z[:,:-(expt.nbins-1)]=x
    # Shift posterior values to end
    z[:,-1]=z[:,ncols-1] # Copy...
    #print z[:,-1]
    #z[z[:,-1].argsort()] # Sort in-place by posterior value
    #z68=
    #print z[:,-1]
    #p68=numpy.percentile(z[:,-1],100-68)
    #p95=numpy.percentile(z[:,-1],100-95)
    #print p68
    #print p95
    #zz=1.0*z
    #print zz.shape
    #print zz[:,-1]

#    sys.exit(0)
    z[:,ncols-1]=0.0     # ...and blank

    # Fetch best-fit parameters and calculate best-fit line
    ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join(outdir,outstem))
    drawml=ana.get_best_fit()['parameters']

    summf=os.path.join(outdir,'1-summary.txt')
    summary=numpy.genfromtxt(summf)[-1,:]
    drawmap=summary[-(ncols+1):-2]

    if False:
        print '--> Calculating *ML* reconstruction'
        drawmap=drawml

    # Convert drawmap into correct units etc.
    power=2.5
    #ymap=expt.evaluate(expt.convertPosterior(drawmap,power))
    ymap=expt.evaluate(drawmap)
    print 'yy',['%e'%y for y in ymap]
    #sys.exit(0)
    if expt.kind!='ppl' or True:
        ymap*=numpy.power(expt.binsMedian/1.0e6,power)
    print 'yy',ymap
    #sys.exit(0)

    for isamp in xrange(nsamp):
        #z[isamp,ncols-1:]=expt.evaluate(expt.convertPosterior(z[isamp,:],power))
        z[isamp,ncols-1:]=expt.evaluate(z[isamp,:])
        if expt.kind!='ppl' or True:
            z[isamp,ncols-1:]*=numpy.power(expt.binsMedian/1.0e6,power)

    # Blanking, 0.0 -> NaN
    z[numpy.where(z==0.0)]='NaN'

    # Save the raw reconstructions
    reconf='recon_raw.txt'
    reconf=os.path.join(outdir,reconf)
    recons=z[:,ncols-1:]
    numpy.savetxt(reconf,recons)

    # Record the 68 and 95 per cent samples
    #print zz[:,-1]
    #zz68=zz[numpy.where(zz[:,-1]>p68)]
    #print zz68.shape
    #zz95=zz[numpy.where(zz[:,-1]>p95)]
    #numpy.savetxt('zz68.txt',zz68)
    #numpy.savetxt('zz95.txt',zz95)
    #from contour_plot import findconfidence
    #zzz=numpy.cumsum(numpy.exp(zz[:,-1]))/2.15430044
    #print zzz
    #a,b,c=findconfidence(zzz)
    #print a,b,c
    #sys.exit(0)

    # Generate stats here...
    s=numpy.zeros((expt.nbins,6))
    s[:,0]=expt.binsMedian

    print '# ibin flux fit low high dlower dupper skew kurtosis'
    for ibin in xrange(expt.nbins):
        x = recons[:,ibin]
        # Remove NaNs from stats vectors
        # http://stackoverflow.com/questions/11620914/removing-nan-values-from-an-array
        x = x[~numpy.isnan(x)]
        #ss=stats.bayes_mvs(x,alpha=0.68)[0]
        #x*=numpy.power(s[ibin,0]/1.0e6,2.5)
        # Select statistic on which to centre reconstruction:
        value_central=ymap[ibin] # MAP
        #value_central=numpy.median(x) # median
        #value_central=peak_confidence(x,bins=10) # peak/mode
        try:
            ss=numpy.zeros(3)
            #ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence(x,alpha=0.68,ret_all=True)
            ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence2(x,alpha=0.95,ret_all=True,\
                                                               value_central=value_central,\
                                                               truncate_edges=True)
        except:
            ss=numpy.nan*numpy.ones(3)

#        tt=peak_confidence(x,bins=10)
        #ss*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print ss[0],tt
#        s[ibin,1]=ss[0]  # median
#        s[ibin,1]=tt     # peak
#        s[ibin,1]=ymap[ibin] # MAP
        s[ibin,1]=value_central
        #print ymap
        s[ibin,2]=value_central-ss[1]  # lower
        s[ibin,3]=ss[2]-value_central  # upper
        s[ibin,4]=stats.skew(x) # skewness
        s[ibin,5]=stats.kurtosis(x) # kurtosis
        print ibin,s[ibin,0],s[ibin,1],dlow,dhigh,ss[1],ss[2],s[ibin,4],s[ibin,5]#,stats.skewtest(x)

    # ...and output to file
    rstatsf='recon_stats.txt'
    rstatsf=os.path.join(outdir,rstatsf)
    hdr='# median_flux_uJy dnds_2p5_Jy1p5srm1 delta_dnds_2p5_lower_Jy1p5srm1 delta_dnds_2p5_upper_Jy1p5srm1 skewness kurtosis'
    fid = open(rstatsf,'w')
    print hdr
    print s

    print '-> Writing stats (i.e. dnds2p5+/- v. S/uJy) to %s' % rstatsf
    fid.write('%s\n'%hdr)
    numpy.savetxt(fid,s)
    fid.close()
    print 'Finished.'

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
