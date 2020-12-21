#!/usr/bin/env python

"""
This is reconstruct_lf.py
June 20 
Eliab based on Jon's reconstruct.py

Usage:
Works for LFs (Sch,DPL,DDPL)
./reconstruct.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
from scipy import stats
import pylab as plt
import pymultinest
from bayestackClasses import countModel
from utils import sqDeg2sr,calculate_confidence,peak_confidence,fetchStats,calculate_confidence2
from lumfuncUtils import get_sbins,get_z,get_dl

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
    #bins=[12.679592, 31.84969509, 80.00281696, 200.9579904, 504.78364938, 1010, 3267.95919972]# 0.4 bins
    bins=numpy.logspace(0,3.01,10)
    bins = list(bins)
    bins.append(2e3)
    bins1 = numpy.arange(17.2,28.2,0.1)
    print 'run expt'
    nbins = len(bins)

    z_m = redshifts[0]
    dl = get_dl(z_m) 
    sbin1 = get_sbins(10**bins1,z_m,dl)*1e6
    expt=countModel(modelFamily,nlaws,settingsf,[dataset],floatNoise,doRedshiftSlices=True,mybins=sbin1)

    #expt=countModel(modelFamily,nlaws,settingsf,dataset,floatNoise,\
     #               doRedshiftSlices=True)#numpy.logspace(0,3.0,20))
                    
    #print 'these are my beens', expt.bins
    #sys.exit()

    f='%spost_equal_weights.dat' % outstem
    #f='%sev.dat'% outstem
    f=os.path.join(outdir,f)
	 # Fetch best-fit parameters and calculate best-fit line
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    
    # Load equally-weighted posterior samples
    x=numpy.genfromtxt(f)
    nsamp=x.shape[0]
    ncols=x.shape[1] # The fifth [seventh] column is the posterior value
    # There must be a better way, but:
    #ncols = 14
    z=numpy.zeros((nsamp,ncols-1+expt.nbins))
    z[:,:-(expt.nbins-1)]=x
    # Shift posterior values to end
    z[:,-1]=z[:,ncols-1] # Copy...
    z[:,ncols-1]=0.0     # ...and blank

    # Fetch best-fit parameters and calculate best-fit line
    ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join(outdir,outstem))
    drawmap=ana.get_best_fit()['parameters']
    #print 'These is MAP',drawml
    #sys.exit()

    #summf=os.path.join(outdir,'1-summary.txt')
    #summary=numpy.genfromtxt(summf)[-1,:]
    #drawmap=summary[-(ncols+1):-2] 
    print 'Now this is MAP',drawmap
    #drawmap = [  1.76792100e+02,   5.24202324e+18,   7.67284256e+24,   2.52157587e+1, 9.63136831e+22,   3.30233162e+00]

    
    

    if True:
        print '--> Calculating *ML* reconstruction'
        drawmap=drawmap

    # Convert drawmap into correct units etc.
    power=2.5
    #drawmap[6]-=0.6
    #ymap=expt.evaluate(expt.convertPosterior(drawmap,power))
    ymap=expt.evaluate(drawmap)[0]#,expt.binsMedian)
    #ymap=expt.realise(drawmap)
    #print 'ymap'
    #print ymap
    #print 'yy',['%e'%y for y in ymap]
    #sys.exit(0)
    
    for isamp in xrange(nsamp):
        #z[isamp,ncols-1:]=expt.evaluate(expt.convertPosterior(z[isamp,:],power))
        #z[isamp,:][6]-=0.6
        z[isamp,ncols-1:]=expt.evaluate(z[isamp,:])[0]#,expt.binsMedian)
        #print z[isamp,:]
        #sys.exit()

    # Blanking, 0.0 -> NaN
    z[numpy.where(z==0.0)]='NaN'
    #sys.exit()

    # Save the raw reconstructions
    reconf='recon_raw.txt'
    reconf=os.path.join(outdir,reconf)
    recons=z[:,ncols-1:]
    #numpy.savetxt(reconf,recons)

    # Generate stats here...
    s=numpy.zeros((len(expt.binsMedian),8))
    s[:,0]=expt.binsMedian

    #print '# ibin flux fit low high dlower dupper skew kurtosis'
    #print '%6s %12s %12s %12s %12s %12s %12s %12s %12s'%('bin', 'lower phi_rms','lower phi', 'l_err','ymap','u_err', 'upper phi','u_phi_rms', '2*rms')
    for ibin in xrange(len(expt.binsMedian)):
        x = recons[:,ibin]
        # Remove NaNs from stats vectors
        # http://stackoverflow.com/questions/11620914/removing-nan-values-from-an-array
        x = x[~numpy.isnan(x)]
        #ss=stats.bayes_mvs(x,alpha=0.68)[0]
        #x*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print bins1[ibin]
        #plt.hist(numpy.log10(x),bins=15)
        #plt.xlabel(r'$\rm{log_{10}[\rho_m(Mpc^{-3} mag^{-1})]}$',fontsize = 20)
        #plt.show()
        #tt=peak_confidence(x,bins=10)

        try:
            ss=numpy.zeros(3)
            ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence(x,alpha=0.95,ret_all=True)
            #ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence2(x,alpha=0.68,ret_all=True,\
            #                                                   value_central=ymap[ibin],\
            #                                                   truncate_edges=True)
        except:
            ss=numpy.nan*numpy.ones(3)
            print "didn't work for ", x
            continue
        #sys.exit()
        #tt=peak_confidence(x,bins=10)
        #ss*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print ss[0],tt
#        s[ibin,1]=ss[0]  # median
#        s[ibin,1]=tt     # peak
        s[ibin,1]=ymap[ibin] # MAP
        #print ymap
        s[ibin,2]=ss[1]  # lower
        s[ibin,3]=ss[2]  # upper
        s[ibin,4]= 2*numpy.std(x)# (s[ibin,3]- s[ibin,2])/(2*numpy.std(x))# skewness
        s[ibin,5]= ss[0]
        
        #print '%6.1f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f'%(bins1[ibin], numpy.log10((ymap[ibin] - s[ibin,4])), numpy.log10(ss[1]),numpy.log10(s[ibin,2]),numpy.log10(ymap[ibin]),numpy.log10(s[ibin,3]) ,numpy.log10(ss[2]), numpy.log10(s[ibin,4]+ ymap[ibin]), numpy.log10(s[ibin,4]) )
        #print ibin,s[ibin,0],s[ibin,1],dlow,dhigh,ss[1],ss[2],s[ibin,4],s[ibin,5]#,stats.skewtest(x)
        #sys.exit()

    # ...and output to file
    #sys.exit()
    rstatsf='recon_stats.txt'
    rstatsf=os.path.join(outdir,rstatsf)
    hdr='# median_flux_uJy dnds_2p5_Jy1p5srm1 delta_dnds_2p5_lower_Jy1p5srm1 delta_dnds_2p5_upper_Jy1p5srm1 skewness kurtosis'
    fid = open(rstatsf,'w')
    #print hdr
    #print s

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
