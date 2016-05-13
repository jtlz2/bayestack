#!/usr/bin/env python

"""
This is plot_counts.py
Jonathan Zwart
December 2015

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
import countUtils
from utils import sqDeg2sr,fetchStats,medianArray
from math import pi

param_file=sys.argv[-1].rstrip('\\')
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

    stokes='S'
    if doPoln: stokes='P'

    expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise,\
                    doPoln=doPoln)
    print expt.data
    #print expt.bins
    #sys.exit(0)
    # Fetch best-fit parameters and calculate best-fit line
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    stats=fetchStats(outdir,expt.parameters,plotTruth)
    SMIN_MAP=stats['S0'][0]
    SMIN_MAP_UPPER=stats['S0'][-2]
    SMIN_MAP_LOWER=stats['S0'][-1]

    # Calculate dn/ds for the data
    j,dn_by_ds_eucl,dn_by_ds_errs,dnbdsb,j=expt.dn_by_ds(return_all=True,\
                                                         data=expt.data)
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
    plt.xlim(1.0e-3,1.0e3)
    #plt.ylim(1.0e-6,0.01)
    #plt.ylim(1.0e11,1.0e13)
    if not doPoln:
        plt.ylim(1.0e-5,100.0)
    else:
        plt.ylim(1.0e-5,10.0)
        #plt.ylim(1.0e-6,100.0)
        #pass

    # Overlay literature values
    if 'en1jvla' in dataset:
        litfile='/home/jtlz2/Dropbox/elaisn1/dezotti/4d8GHz.dat'
        lit=numpy.genfromtxt(litfile)
        xlit=lit[:,0]; ylit=lit[:,1]
        dylit_up=lit[:,2]; dylit_down=lit[:,3]
        if not doPoln:
            plt.errorbar(10**(6+xlit),ylit,yerr=[dylit_up,dylit_down],\
                         fmt='k.',alpha=0.25,label='literature 4.8 GHz')

    # Overlay Massardi models
    # These are reported as being nS^0, but in fact they are nS^2.5 :
    powerfix=0.0
    mmfile='/home/jtlz2/Dropbox/elaisn1/dezotti/model_all_4d8GHz.dat'
    mm=numpy.genfromtxt(mmfile)
    xmm=mm[:,0]
    ymm=10**(2.5*powerfix*xmm+mm[:,1]); y2mm=10**(2.5*powerfix*xmm+mm[:,3])
    if not doPoln and 'en1jvla' in dataset:
        plt.plot(10**(6+xmm),ymm,'b-',alpha=0.25,label='Massardi 4.8 GHz (radio)')
        plt.plot(10**(6+xmm),y2mm,'g-',alpha=0.25,label='Massardi 4.8 GHz (sbn)')
        plt.plot(10**(6+xmm),ymm+y2mm,'r-',alpha=0.25,\
                 label='Massardi 4.8 GHz (radio+sbn)')

    # Polarization counts
    if doPoln and 'en1jvla' in dataset:
        plt.text(0.1e2,2.0,'4.8 GHz poln',color='k',fontsize=24)
        alpha_r=-0.7
        r=(4.8/1.4)**alpha_r
        # Overlay Stil NVSS 1.4-GHz counts
        # !!Need to do a conversion to 5 GHz!!
        nvssf='/home/jtlz2/Dropbox/elaisn1/stil/nvss_stil.txt'
        nvss=numpy.genfromtxt(nvssf)
        xnvss=1.0e3*nvss[:,0]; ynvss=r*nvss[:,1]
        plt.plot(xnvss,ynvss,'g+',label='NVSS (Stil+ 2014)')

        # Overlay Banfield ELAIS-N1 1.4-GHz counts
        pass

        # Overlay Hales ATLAS-DR2 1.4-GHz counts
        halesf1='/home/jtlz2/Dropbox/elaisn1/hales/hales_A5.txt'
        hales1=numpy.genfromtxt(halesf1)
        xhales1=1.0e3*hales1[:,2]; yhales1=r*hales1[:,15]
        dyhales1_up=r*hales1[:,16]; dyhales1_down=r*hales1[:,17]
        plt.errorbar(xhales1,yhales1,yerr=[-dyhales1_down,dyhales1_up],\
                     fmt='g.',label='ATLAS / CDF-S (Hales+ 2014)')

        halesf2='/home/jtlz2/Dropbox/elaisn1/hales/hales_A6.txt'
        hales2=numpy.genfromtxt(halesf2)
        xhales2=1.0e3*hales2[:,2]; yhales2=r*hales2[:,15]
        dyhales2_up=r*hales2[:,16]; dyhales2_down=r*hales2[:,17]
        plt.errorbar(xhales2,yhales2,yerr=[-dyhales2_down,dyhales2_up],\
                     fmt='m.',label='ATLAS / ELAIS-S1 (Hales+ 2014)')

        # Plot SKADS 5-GHz poln counts
        skadspoldir='/home/jtlz2/Dropbox/elaisn1/skadspol'
        skadsf={0:'All_pcounts_5GHz.dat',1:'FRI_pcounts_5GHz.dat',2:'FRII_pcounts_5GHz.dat',\
                3:'RQQ_pcounts_5GHz.dat',4:'NG_pcounts_5GHz.dat'}
        styles={0:'r-',1:'c--',2:'m--',3:'y--',4:'b--'}
        skads={}
        for sfnum in skadsf.keys():
            sf=os.path.join(skadspoldir,skadsf[sfnum])
            skads[sfnum]=numpy.genfromtxt(sf)
            xskads=1.0e9*skads[sfnum][:,0];
            nonzero=numpy.where(xskads!=0.0)
            yskads=numpy.power(skads[sfnum][:,0],2.5)*skads[sfnum][:,2]
            label=skadsf[sfnum].split('_')[0]
            plt.plot(1.0e-3*xskads[nonzero],yskads[nonzero],styles[sfnum],label='SKADS (%s)'%label)

    if not doPoln and 'goodsnx' in dataset:
            # Eric's total-intensity 5-sig source-detection catalogue
            fivesigcat='/home/jtlz2/Dropbox/bayestack/eric/gncat_fluxes_5sig.txt'
            fluxes=numpy.genfromtxt(fivesigcat)
            bins=numpy.array([3.5,6.0,8.0,10.0,15.0,30.0,\
                              50.0,75.0,100.0,370.0])
            counts=numpy.histogram(fluxes,bins=bins)[0]
            area5sig=pi*(6.5/60.0)**2
            # Calculate differential counts
            idl_style=False
            dn_by_ds,dn_by_ds_eucl,dn_by_ds_errs,dn_by_ds_b,dn_by_ds_b_errs=\
                countUtils.calculateDnByDs(1.0e-6*bins,counts,\
                                           idl_style=idl_style,return_all=True)
            plt.errorbar(medianArray(bins),dn_by_ds_eucl/area5sig/sqDeg2sr,\
                         fmt='k+',yerr=dn_by_ds_errs/area5sig/sqDeg2sr,label='detected')
            print bins
            print counts
    # Overlay > 5 sigma values
    # See http://localhost:8888/notebooks/Dropbox/bayestack/notebooks/deepcounts.ipynb
    if not doPoln and 'en1jvla' in dataset:
        plt.text(0.1e2,20.0,'4.8 GHz',color='k',fontsize=24)

        fivesigf='/home/jtlz2/Dropbox/elaisn1/jvla_5GHz_counts.txt'
        fivesig=numpy.genfromtxt(fivesigf)
        #print fivesig[:,1]
        cutoff2=numpy.where(fivesig[:,2]>5.0) # muJy
        plt.errorbar(fivesig[:,2][cutoff2],fivesig[:,5][cutoff2],fmt='b*',fillstyle='none',\
                 markerfacecolor='w',yerr=fivesig[:,6][cutoff2],label='uncorrected')
        plt.errorbar(fivesig[:,2][cutoff2],fivesig[:,3][cutoff2],fmt='b*',\
                     yerr=fivesig[:,4][cutoff2],label='corrected')
        cutoff=numpy.where(fivesig[:,2]>1.0e-2) # muJy
        plt.errorbar(fivesig[:,2][cutoff],fivesig[:,7][cutoff],fmt='r.-',\
                 yerr=fivesig[:,8][cutoff],label='SKADS')
    elif 'en1jvla' in dataset:
        # Plot direct directions
        directf='/home/jtlz2/Dropbox/elaisn1/russ/JVLA_5GHz_direct_detns.txt'
        direct=numpy.genfromtxt(directf)
        xdirect=1.0e3*direct[:,2]; ydirect=direct[:,5]; yerrdirect=direct[:,6]
        plt.errorbar(xdirect,ydirect,yerr=yerrdirect,fmt='k.',label='direct detections')
    elif 'goodsnx' in dataset:
        plt.text(0.1e2,20.0,'10.0 GHz',color='k',fontsize=24)

    # Plot the data
    #print expt.binsMedian
    #print dn_by_ds_eucl
    #print SURVEY_AREA
    #print expt.survey.SURVEY_AREA
    #print expt.data
    #sys.exit(0)
#    plt.errorbar(expt.binsMedian[numpy.where(dn_by_ds_eucl>0)],\
#                 dn_by_ds_eucl[numpy.where(dn_by_ds_eucl>0)]/SURVEY_AREA/sqDeg2sr,\
#                 fmt='b+',yerr=dn_by_ds_errs[numpy.where(dn_by_ds_eucl>0)]\
#                 /SURVEY_AREA/sqDeg2sr,label='data')

    # Plot the reconstruction
    xrecon=s[:-1,0]; yrecon=s[:-1,1]
    yrecon_down=s[:-1,2]; yrecon_up=s[:-1,3]

    smask = (xrecon>SMIN_MAP)
    plt.fill_between(xrecon[smask],(yrecon[smask]-yrecon_down[smask]),\
                     (yrecon[smask]+yrecon_up[smask]),color='k',alpha=0.2)
    plt.errorbar(xrecon[smask],yrecon[smask],fmt='k',label='MAP estimate')
    plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
    plt.axvline(SMIN_MAP,color='r',alpha=1.0)
    #plt.axvline(SMIN_MAP_UPPER,color='r',alpha=0.2)
    #plt.axvline(SMIN_MAP_LOWER,color='r',alpha=0.2)

    plt.axvspan(SMIN_MAP_LOWER,SMIN_MAP_UPPER,alpha=0.1,color='r',\
                label=r'$68\%\,\mathrm{interval}\,'+stokes+'_{min}$'+' (%3.1f nJy)'%(1.0e3*SMIN_MAP))
    #plt.axvline(1.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)
    #plt.axvline(5.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),color='g',alpha=0.2)
    if True:
        plt.axvspan(1.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),\
                5.0*SURVEY_NOISE/numpy.sqrt(expt.nsrc),alpha=0.1,color='g',\
                label=r'$\left(1-5\right)\sigma/\sqrt{N=%s}$'%expt.nsrc)

    #plt.text(SURVEY_NOISE,0.16,'1 sigma',rotation=90,color='b',alpha=0.5)
    #plt.text(5.0*SURVEY_NOISE,0.16,'5 sigma',rotation=90,color='b',alpha=0.5)

    f=10
    legend=plt.legend(loc='upper left',prop={'size':9},frameon=True,numpoints=1)
    frame = legend.get_frame()
    frame.set_facecolor('white')

    plt.xlabel('$%s/\mu\mathrm{Jy}$'%stokes,fontsize=f)
    plt.ylabel('$n%s^{2.5}/\mathrm{sr}^{-1}\mathrm{Jy}^{1.5}$'%stokes,fontsize=f)
    plt.tight_layout()

    plotf='%s/recon_plot_%s.pdf' % (outdir,run_num)
    plt.savefig(plotf)
    print '-> Run: open %s' % plotf
    plt.close()

    # Now make the linear plots

    #fig = plt.figure()
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.xlim(1.0,2000.0)
    #plt.xlim(1.0e-3,1.0e3)
    #plt.ylim(1.0e-6,0.01)
    #plt.ylim(1.0e11,1.0e13)
    #if not doPoln:
    #    plt.ylim(1.0e-5,100.0)
    #else:
    #    plt.ylim(1.0e-5,10.0)
    #    #plt.ylim(1.0e-6,100.0)
    #    #pass

    #plotf2='%s/recon_linear_%s.pdf' % (outdir,run_num)
    #plt.savefig(plotf2)
    #print '-> Run: open %s' % plotf2
    #plt.close()
    
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
