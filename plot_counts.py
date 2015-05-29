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
from utils import *
import contour_plot
from bayestackClasses import countModel

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

    s=pylab.loadtxt('%s/recon_stats.txt'%outdir)

    fig = plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(bins[0],20.0*SURVEY_NOISE)
    #plt.ylim(0.5,1.0e3)
    plt.xlabel('S/$\mu$Jy')
    plt.ylabel('Number of objects')

    xrecon=s[:,0]; yrecon=s[:,1]
    yrecon_down=s[:,2]; yrecon_up=s[:,3]
    plt.fill_between(xrecon,yrecon-yrecon_down,yrecon+yrecon_up,color='k',alpha=0.2)
    plt.errorbar(xrecon,yrecon,fmt='k',label='MAP estimate')
    plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
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
