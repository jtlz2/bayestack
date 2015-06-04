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
from utils import sqDeg2sr,calculate_confidence,peak_confidence,medianArray

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
    expt=countModel(modelFamily,nlaws,settingsf,dataset,floatNoise)

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
    ymap=expt.evaluate(expt.convertPosterior(drawmap,power))
    if expt.kind=='poly':ymap*=numpy.power(expt.binsMedian/1.0e6,2.5)
    #print drawmap
    #print ymap
    #sys.exit(0)
    for isamp in xrange(nsamp):
        z[isamp,ncols-1:]=expt.evaluate(expt.convertPosterior(z[isamp,:],power))
        if expt.kind=='poly':z[isamp,ncols-1:]*=numpy.power(expt.binsMedian/1.0e6,2.5)

    # Blanking, 0.0 -> NaN
    z[numpy.where(z==0.0)]='NaN'

    # Save the raw reconstructions
    reconf='recon_raw.txt'
    reconf=os.path.join(outdir,reconf)
    recons=z[:,ncols-1:]
    numpy.savetxt(reconf,recons)

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

        try:
            ss=numpy.zeros(3)
            ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence(x,alpha=0.68,ret_all=True)
        except:
            ss=numpy.nan*numpy.ones(3)
        tt=peak_confidence(x,bins=10)
        #ss*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print ss[0],tt
#        s[ibin,1]=ss[0]  # median
        #s[ibin,1]=tt     # peak
        s[ibin,1]=ymap[ibin] # MAP
        #print ymap
        s[ibin,2]=ss[0]-ss[1]  # lower
        s[ibin,3]=ss[2]-ss[0]  # upper
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
