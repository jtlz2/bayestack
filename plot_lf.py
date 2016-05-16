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
from countUtils import calculateDnByDs,medianArray
from dnds_lumfunc import get_Lbins,get_sbins,get_dl
from plat import open_anytype, calc_Vmax
from matplotlib.ticker import AutoMinorLocator
#import sourcecounts1

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file

@profile
def main():

    """
    """

    print 'Settings file is %s' % setf

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
    recon_expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise,\
        mybins=numpy.logspace(0.0,3.0,200))
    print expt.data
    print expt.bins
    # Fetch best-fit parameters and calculate best-fit line
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    stats=fetchStats(outdir,expt.parameters,plotTruth)
    stats=fetchStats(outdir,expt.parameters,plotTruth)
    SMIN_MAP=stats['LMIN'][0]
    SMIN_MAP_UPPER= stats['LMIN'][-2]
    SMIN_MAP_LOWER=stats['LMIN'][-1]
    #print stats
    #sys.exit()

    
    #SMIN_MAP=stats['S0'][0]
    #SMIN_MAP_UPPER= stats['S0'][-2]
    #SMIN_MAP_LOWER=stats['S0'][-1]

    # Calculate dn/ds for the data
    j,dn_by_ds_eucl,dn_by_ds_errs,dnbdsb,j=expt.dn_by_ds(return_all=True,data=expt.data)
    if dataset=='sdss':
        d=numpy.genfromtxt('%s/%s'%(outdir,datafile[4:]))
    elif 'sim' in dataset:
        d=numpy.genfromtxt(os.path.join(outdir,'sim.txt'))

    # Load the reconstruction
    s=pylab.loadtxt('%s/recon_stats.txt'%outdir)
    #z = [0, 0.5, 1 , 1.5, 2, 2.3, 2.6, 3, 3.5, 4] #old
    #z = [ 0.7, 1 , 1.35, 1.7, 2, 2.3, 2.6, 3, 3.5, 4]# older
    z =[0.2, 0.45, 0.7, 1.0, 1.3, 1.6, 1.85, 2.15, 2.35, 2.55, 2.85, 3.15, 3.5]#oldest ;)
    num = int(datafile[4:][-5])
    print num
    
    
    z_min= z[int(num) -1]
    z_max= z[int(num)]
    z_m = (z_min + z_max)/2
    print z_min,z_max
    xrecon=s[:-1,0]; yrecon=s[:-1,1]
    yrecon_down=s[:-1,2]; yrecon_up=s[:-1,3]
    dl = get_dl(z_m)
    SMIN = get_sbins(SMIN_MAP,z_m,dl)*1e6
    Lbins=get_Lbins(xrecon[numpy.where(xrecon>SMIN)],z_m,dl)
    phi = yrecon[numpy.where(xrecon>SMIN)]
    lin_phi = numpy.log10(phi)
    Lbins = numpy.log10(Lbins)
    
    
    zef,l =open_anytype('../First_dr7d_1_8arc.txt',(2,9),[z_min,z_max]) # dr7
    bins = numpy.arange(23.6,29.2,0.4)
    rho_1,er1 = calc_Vmax(zef,l,bins,z_max,z_min,SURVEY_AREA*sqDeg2sr)
    
    fig,ax = plt.subplots()
    plt.errorbar(Lbins,lin_phi,fmt='x',label='MAP')
    plt.errorbar(bins,rho_1,yerr=er1,fmt = 'ob',label='Detected')
    plt.legend()
    
    plt.xlabel('1.4 GHz log[L(W/Hz)]',fontsize = 18)
    plt.ylabel(r'log[$\rho_m$(mpc$^{-3}$mag$^{-1}$)]',fontsize = 18)
    plt.tick_params(axis='both',which = 'major', labelsize=15,width =2)
    plt.tick_params(axis='both',which = 'minor', labelsize=12, width=1)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())

    
    plt.show()
    
	
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
