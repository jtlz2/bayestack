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
    SMIN_MAP_UPPER= stats['S0'][-2]
    SMIN_MAP_LOWER=stats['S0'][-1]

    # Calculate dn/ds for the data
    j,dn_by_ds_eucl,dn_by_ds_errs,dnbdsb,j=expt.dn_by_ds(return_all=True,data=expt.data)
    if dataset=='sdss':
        d=numpy.genfromtxt('%s/%s'%(outdir,datafile[4:]))
    elif 'sim' in dataset:
        d=numpy.genfromtxt(os.path.join(outdir,'sim.txt'))

    #for b in range(expt.nbins):
    #    print b,expt.binsMedian[b],dnbdsb[b]
    #sys.exit(0)
    # Load the reconstruction
    s=pylab.loadtxt('%s/recon_stats.txt'%outdir)
    #z = [0, 0.5, 1 , 1.5, 2, 2.3, 2.6, 3, 3.5, 4] #old
    #z = [ 0.7, 1 , 1.35, 1.7, 2, 2.3, 2.6, 3, 3.5, 4]# older
    z =[0.2, 0.45, 0.7, 1.0, 1.3, 1.6, 1.85, 2.15, 2.35, 2.55, 2.85, 3.15, 3.5]#oldest ;)
    models = ['A','B','C','D',"A'","B'","C'","D'"]
    #models = ["A'","B'","C'","D'",'A','B','C','D']
    a = ['a','i','q']
    b = ['b','j','r']
    c = ['c', 'k','s']
    d = ['d','l','t']
    e = ['e','m','u']
    f = ['f','n','v']
    g = ['g','o','w']
    h = ['h','p','x']
    letters =[a,b,c,d,e,f,g,h]
    z_max,z_min = 20.,0.
    z_min= z[int(num) -1]
    z_max= z[int(num)]
    print z_min,z_max
    	 
    #z_max,z_min =2.5,1.8
    scube=pylab.loadtxt('../skads_cat2.el')
    area_skads = 3.365883939231586
    
    f2 = numpy.loadtxt('../catalogue_dr12d_1_8.txt') # dr12_1d_18.el
    flux2 = f2[numpy.where((f2[:,4]> z_min) & (f2[:,4]< z_max) )][:,5]*1e-3
    
    #f2 = numpy.loadtxt('../dr12_1d_18.el')
    #flux2 = f2[:,-7]*1e-3# Condon sample
	
    print len(scube)
    scube = scube[numpy.where((scube[:,2]> z_min) & (scube[:,2]< z_max) )] #getting the write sample    
    rqq = scube[:,-3][numpy.where(scube[:,-1] == 1)]
    print numpy.shape(rqq)	
    rlq = scube[:,-3][numpy.where(numpy.logical_or(scube[:,-1] == 2,scube[:,-1] ==3)) ]
    
    h1 = numpy.histogram(scube[:,-3],bins=numpy.logspace(-6,1,20))
    print numpy.shape(scube)
    x1    = medianArray(h1[1]*1e6)
    y1    = calculateDnByDs(h1[1],h1[0],eucl=True,idl_style=False) /area_skads/sqDeg2sr
    
    h2 = numpy.histogram(rqq,bins=numpy.logspace(-6,1,20))
    x2    = medianArray(h2[1]*1e6)
    y2    = calculateDnByDs(h2[1],h2[0],eucl=True,idl_style=False) /area_skads/sqDeg2sr
    
    h3 = numpy.histogram(rlq,bins=numpy.logspace(-6,1,20))
    x3    = medianArray(h3[1]*1e6)
    y3    = calculateDnByDs(h3[1],h3[0],eucl=True,idl_style=False) /area_skads/sqDeg2sr
    
    h4    = numpy.histogram(flux2,bins=numpy.logspace(-6,1,30))
    x4    = medianArray(h4[1]*1e6)
    y4    = calculateDnByDs(h4[1],h4[0],eucl=True,idl_style=False) /SURVEY_AREA/sqDeg2sr
    err4  = calculateDnByDs(h4[1],h4[0],errors=True,idl_style=False)/SURVEY_AREA/sqDeg2sr
    
    #h0 = numpy.histogram(data,bins=numpy.logspace(-6,1,50))
    print max(rqq)*1e6
    #x0 = medianArray(h0[1]*1e6)
    #y0 = calculateDnByDs(h0[1],h0[0],eucl=True,idl_style=False)/SURVEY_AREA/sqDeg2sr
    #print y0
  

    fig = plt.figure()
    plt.xscale('log')
    plt.yscale('log')
    #plt.xlim(1.0,2000.0)
    #plt.ylim(1.0e-6,0.01)
    #plt.ylim(1.0e11,1.0e13)

    # Plot the data
    #print expt.binsMedian
    #print dn_by_ds_eucl
    #print expt.data
    #plt.plot(x0,y0,'rx',label='Extracted fluxes')
    plt.plot(x1,y1/70,label='Wilman et. al 2008, All sources')
    plt.plot(x2,y2/70,label='Wilman et. al 2008, radio-quiet')
    plt.plot(x3,y3/70,label='Wilman et. al 2008, radio-loud')
    plt.errorbar(x4,y4,err4,fmt='.',label='Detected sources')
    #plt.errorbar(expt.binsMedian[numpy.where(dn_by_ds_eucl>0)],\
    #             dn_by_ds_eucl[numpy.where(dn_by_ds_eucl>0)]/SURVEY_AREA/sqDeg2sr,\
    #             fmt='b+',yerr=dn_by_ds_errs[numpy.where(dn_by_ds_eucl>0)]\
    #             /SURVEY_AREA/sqDeg2sr,label='data')

    # Plot the reconstruction
    xrecon=s[:-1,0]; yrecon=s[:-1,1]
    yrecon_down=s[:-1,2]; yrecon_up=s[:-1,3]

    plt.fill_between(xrecon[xrecon>SMIN_MAP],\
                     (yrecon[xrecon>SMIN_MAP]-yrecon_down[xrecon>SMIN_MAP])\
                     /SURVEY_AREA/sqDeg2sr,(yrecon[xrecon>SMIN_MAP]+\
                     yrecon_up[xrecon>SMIN_MAP])/SURVEY_AREA/sqDeg2sr,\
                     color='k',alpha=0.2)
    #plt.errorbar(xrecon,(yrecon)/SURVEY_AREA/sqDeg2sr,fmt='k',label='MAP estimate')
    plt.errorbar(xrecon[xrecon>SMIN_MAP],(yrecon[xrecon>SMIN_MAP])/SURVEY_AREA/sqDeg2sr,fmt='k',label='MAP estimate')
    plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(SMIN_MAP,color='r',alpha=1.0)
    #plt.axvline(SMIN_MAP_UPPER,color='r',alpha=0.2)
    #plt.axvline(SMIN_MAP_LOWER,color='r',alpha=0.2)
    plt.axvspan(SMIN_MAP_LOWER,SMIN_MAP_UPPER,alpha=0.1,color='r',\
                label='$95\%\,\mathrm{interval}\,S_{min}$')
    #plt.axvline(1.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)
    #plt.axvline(5.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)
    plt.axvspan(1.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),\
                5.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),alpha=0.1,color='g',\
                label=r'$\left(1-5\right)\sigma/\sqrt{N=%s}$'%expt.nsrc)

    #plt.text(SURVEY_NOISE,0.16,'1 sigma',rotation=90,color='b',alpha=0.5)
    #plt.text(5.0*SURVEY_NOISE,0.16,'5 sigma',rotation=90,color='b',alpha=0.5)
    #plt.text(0.15,10**-2,'%.2f < z < %.2f'%(z_min,z_max),alpha=0.5 )
    
    for i in range(len(models)):
    	if param_file[-1] in letters[i]:
    		plt.text(40,10,models[i],fontsize=18)
    plt.text(5e4,10**1,'%.2f < z < %.2f'%(z_min,z_max),alpha=0.5 )
    f=16
    legend=plt.legend(loc='lower right',prop={'size':11},frameon=False,numpoints=1)
    #frame = legend.get_frame()
    #frame.set_facecolor('white')
    plt.xlim(1,1e7)
    plt.ylim(1e-6,1e2)
    plt.xlabel('$S/\mu\mathrm{Jy}$',fontsize=f)
    plt.ylabel('$nS^{2.5}/\mathrm{sr}^{-1}\mathrm{Jy}^{1.5}$',fontsize=f)

    plotf='%s/recon_plot_%s.pdf' % (outdir,run_num)
    #plt.show()
    plt.savefig(plotf)
    print '-> Look in %s' % plotf 

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
