#!/usr/bin/env python

"""
This is plot_counts.py
Jonathan Zwart
May 2014

Usage:

./plot_counts.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
if os.getenv('PBS_O_HOST') not in ['baltasar']:
    from matplotlib import pyplot as plt
import pylab
from profile_support import profile
import pymultinest
from bayestackClasses import countModel
from utils import sqDeg2sr,fetchStats

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file

#-------------------------------------------------------------------------------


@profile
def main():

    """
    """

    print 'Settings file is %s' % setf

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
    print expt.data
    print expt.bins
    # Fetch best-fit parameters and calculate best-fit line
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    stats=fetchStats(outdir,expt.parameters,plotTruth)
    SMIN_MAP=stats['S0'][0]
    SMIN_MAP_UPPER=stats['S0'][-2]
    SMIN_MAP_LOWER=stats['S0'][-1]

    # Calculate dn/ds for the data
    j,dn_by_ds_eucl,dn_by_ds_errs,dnbdsb,j=expt.dn_by_ds(return_all=True,data=expt.data)
    if dataset=='sdss':
        d=numpy.genfromtxt('%s/sdss_dr12s1.txt'%outdir)
    elif 'sim' in dataset:
        d=numpy.genfromtxt(os.path.join(outdir,'sim.txt'))

    #for b in range(expt.nbins):
    #    print b,expt.binsMedian[b],dnbdsb[b]
    #sys.exit(0)
    # Load the reconstruction
    s=pylab.loadtxt('%s/recon_stats.txt'%outdir)

    fig = plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(1.0,2000.0)
    #plt.ylim(1.0e-6,0.01)
    #plt.ylim(1.0e11,1.0e13)

    # Plot the data
    #print expt.binsMedian
    #print dn_by_ds_eucl
    #print SURVEY_AREA
    #print expt.data
    plt.errorbar(expt.binsMedian[numpy.where(dn_by_ds_eucl>0)],\
                 dn_by_ds_eucl[numpy.where(dn_by_ds_eucl>0)]/SURVEY_AREA/sqDeg2sr,\
                 fmt='b+',yerr=dn_by_ds_errs[numpy.where(dn_by_ds_eucl>0)]\
                 /SURVEY_AREA/sqDeg2sr,label='data')

    # Plot the reconstruction
    xrecon=s[:-1,0]; yrecon=s[:-1,1]
    yrecon_down=s[:-1,2]; yrecon_up=s[:-1,3]

    plt.fill_between(xrecon[xrecon>SMIN_MAP],\
                     (yrecon[xrecon>SMIN_MAP]-yrecon_down[xrecon>SMIN_MAP]),\
                     (yrecon[xrecon>SMIN_MAP]+yrecon_up[xrecon>SMIN_MAP]),\
                     color='k',alpha=0.2)
    plt.errorbar(xrecon[xrecon>SMIN_MAP],yrecon[xrecon>SMIN_MAP],fmt='k',label='MAP estimate')
    plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(SMIN_MAP,color='r',alpha=1.0)
    #plt.axvline(SMIN_MAP_UPPER,color='r',alpha=0.2)
    #plt.axvline(SMIN_MAP_LOWER,color='r',alpha=0.2)
    plt.axvspan(SMIN_MAP_LOWER,SMIN_MAP_UPPER,alpha=0.1,color='r',\
                label='$68\%\,\mathrm{interval}\,S_{min}$')
    #plt.axvline(1.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)
    #plt.axvline(5.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)
    plt.axvspan(1.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),\
                5.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),alpha=0.1,color='g',\
                label=r'$\left(1-5\right)\sigma/\sqrt{N=%s}$'%expt.nsrc)

    #plt.text(SURVEY_NOISE,0.16,'1 sigma',rotation=90,color='b',alpha=0.5)
    #plt.text(5.0*SURVEY_NOISE,0.16,'5 sigma',rotation=90,color='b',alpha=0.5)

    f=10
    legend=plt.legend(loc='upper left',prop={'size':12},frameon=True,numpoints=1)
    frame = legend.get_frame()
    frame.set_facecolor('white')
    plt.xlabel('$S/\mu\mathrm{Jy}$',fontsize=f)
    plt.ylabel('$nS^{2.5}/\mathrm{sr}^{-1}\mathrm{Jy}^{1.5}$',fontsize=f)

    plotf='%s/recon_plot_%s.pdf' % (outdir,run_num)
    plt.savefig(plotf)
    print '-> Look in %s' % plotf 

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
