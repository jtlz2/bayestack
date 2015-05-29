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

    # Fetch best-fit parameters and calculate best-fit line
    #SMIN_MAP=drawmap[expt.parameters.index('S0')]
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    stats=fetchStats(outdir,expt.parameters,plotTruth)
    SMIN_MAP=stats['S0'][0]
    SMIN_MAP_UPPER=stats['S0'][-2]
    SMIN_MAP_LOWER=stats['S0'][-1]

    s=pylab.loadtxt('%s/recon_stats.txt'%outdir)

    fig = plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(1.0,2000.0)
    plt.ylim(0.01,10.0)

    xrecon=s[:,0]; yrecon=s[:,1]
    yrecon_down=s[:,2]; yrecon_up=s[:,3]
    plt.fill_between(xrecon,(yrecon-yrecon_down)/SURVEY_AREA/sqDeg2sr,\
                     (yrecon+yrecon_up)/SURVEY_AREA/sqDeg2sr,color='k',alpha=0.2)
    plt.errorbar(xrecon,yrecon/SURVEY_AREA/sqDeg2sr,fmt='k',label='MAP estimate')
    plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(SMIN_MAP,color='r',alpha=1.0)
    plt.axvline(SMIN_MAP_UPPER,color='r',alpha=0.2)
    plt.axvline(SMIN_MAP_LOWER,color='r',alpha=0.2)
    plt.axvline(1.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)
    plt.axvline(5.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)

    #plt.text(SURVEY_NOISE,0.16,'1 sigma',rotation=90,color='b',alpha=0.5)
    #plt.text(5.0*SURVEY_NOISE,0.16,'5 sigma',rotation=90,color='b',alpha=0.5)

    f=10
    plt.legend(loc='upper left',prop={'size':12},frameon=False,numpoints=1)
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
