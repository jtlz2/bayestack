#!/usr/bin/env python

"""
This is binner.py
Jonathan Zwart
May 2015

Bin data set for bayestack.py/lumfunc.py

Usage:

./binner.py SETTINGS_FILE.py

"""

import os,sys,shutil
import ctypes
import time,threading
import importlib
import numpy
from math import exp,log,log10,sqrt,isinf,isnan
from scipy import integrate,stats
from profile_support import profile
from utils import *
import lumfunc
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

__name_cached=__name__
param_file=sys.argv[-1]

#setf='%s' % param_file.split('.')[-2]
# Import the settings variables
if __name__=='__main__':
    #setf='binner_settings'
    setf=sys.argv[-1].split('.')[-2]
else:
    setf='%s'%param_file

set_module=importlib.import_module(setf)
globals().update(set_module.__dict__)
__name__=__name_cached

#try:
#    execfile(param_file)
#except IOError:
#    from settings import *


#-------------------------------------------------------------------------------

@profile
def main():
    """
    """

    #global OUTPUT,DUMP
    global CORR_BINS,bins
    
    #CORR_BINS=numpy.ones(len(bins)-1)
    # Set up some fake (unity) correction factors or from frac.py
    #if BIN_CAT_FORM not in [1,3]:
    #CORR_BINS=numpy.ones(len(bins)-1)

    if CORR_BINS is None:
        CORR_BINS=numpy.ones(len(bins)-1)

    print 'Reading from %s' % BIN_CAT

    # COSMOS-VLA data [from Ketron]:
    if BIN_CAT_FORM==0:
        cat=numpy.genfromtxt(BIN_CAT) # mJy in SURVEY_AREA sq. deg.
        cat[:,BIN_COL] *= 1000.0 # mJy -> uJy in SURVEY_AREA sq. deg.
    
        #smap=cat[:,0] *1000.0 # uJy
        #scat=cat[:,1] *1000.0 # uJy

    # VIDEO catalogue form from stacker.py
    elif BIN_CAT_FORM==1 or BIN_CAT_FORM==4:
        cat=numpy.genfromtxt(BIN_CAT) # uJy in SURVEY_AREA sq. deg.
        if len(CORR_BINS) != len(bins)-1:
            print '**Binning corrections mismatch %s' % BIN_CAT,bins,CORR_BINS
            sys.exit(1)

    # VVDF Bondi 2003 catalogue from web
    elif BIN_CAT_FORM==2:
        cat=numpy.genfromtxt(BIN_CAT)
        cat[:,BIN_COL] *= 1000.0
        # 5 sigma:
        #cat=cat[numpy.where((cat[:,BIN_COL]/cat[:,BIN_COL+1])>0.0)]

    # LR-matched Bondi 2003 catalogue from Kim
    elif BIN_CAT_FORM==3:
        cat=numpy.genfromtxt(BIN_CAT)
        cat[:,BIN_COL] *= 1000.0
        #if bins != numpy.array([0.08,0.12,0.18,\
        #                     0.27,0.41,0.61,0.91,1.37,2.05,3.08,4.61,\
        #                     6.92,10.38,15.57,119.23]):
        if len(CORR_BINS) != len(bins)-1:
            print '**Binning corrections mismatch %s' % BIN_CAT,bins,CORR_BINS
            sys.exit(1)    

    # Eliab's SDSS catalogue
    elif BIN_CAT_FORM==6:
        cat=numpy.genfromtxt(BIN_CAT)
        cat[:,BIN_COL] *= 1.0e6
        # 5 sigma:
        #cat=cat[numpy.where((cat[:,BIN_COL]/cat[:,BIN_COL+1])>0.0)]

    # Correct the fluxes for the resolution bias
    if CORR_RESOLUTION is not None:
        print '--> Corrected fluxes for resolution bias (x %f)' % CORR_RESOLUTION
        cat[:,BIN_COL] *= CORR_RESOLUTION

    # Optionally threshold the catalogue
    if BIN_CAT_CLIP is not None:
        cat=cat[numpy.where(cat[:,BIN_COL_CLIP]>BIN_CAT_CLIP)]
        #cat=cat[numpy.where(numpy.abs(cat[:,BIN_COL_CLIP])>BIN_CAT_CLIP)]
    print 'S/uJy: %f -> %f' % (numpy.min(cat[:,BIN_COL]),numpy.max(cat[:,BIN_COL]))
    # Bin the raw data...
    # bins is in uJy; cat is also now in uJy
    counts=numpy.histogram(cat[:,BIN_COL],bins=bins)[0] # uJy in SURVEY_AREA sq. deg.

    # Calculate dn/ds and convert to Jy^-1 sr^-1 and Jy^1.5 sr^-1
    # 6 = -6 x -1
    # 9 = -6 x 1.5
    # Since the counts from the cat are for SURVEY_AREA sq. deg.,
    #           need to correct -> 1 sq. deg. -> 1 sr
    # ACTUALLY - convert this below when writing the table file
    #   -> These are for SURVEY_AREA sq. deg.
    # idl_style or not....? - NO
    idl_s=False
    dn_by_ds=lumfunc.calculateDnByDs(1.0e-6*bins,counts,eucl=False,idl_style=idl_s)
        #/(sqDeg2sr*SURVEY_AREA)
        #*10**6.0/(SURVEY_AREA*sqDeg2sr)
    dn_by_ds_eucl=lumfunc.calculateDnByDs(1.0e-6*bins,counts,eucl=True,idl_style=idl_s)
        #/(sqDeg2sr*SURVEY_AREA)
    #  *10.0**-9.0/(SURVEY_AREA*sqDeg2sr)
    dn_by_ds_errs=lumfunc.calculateDnByDs(1.0e-6*bins,counts,errors=True,idl_style=idl_s)

    # Determine the bin medians, as usual
    median_bins=medianArray(bins) # uJy
    print median_bins


    # ...and write the binned data to file
    if BOUT_CAT is not None:
        outputf=BOUT_CAT
        s=open(outputf,'w')
        header='# bin_low_uJy bin_high_uJy bin_median_uJy Ns_tot_obs dnds_srm1Jym1 dnds_eucl_srm1Jy1p5 delta_dnds_eucl_lower_srm1Jy1p5 delta_dnds_eucl_upper_srm1Jy1p5 corr Ncgts_degm2 dNcgts_lower_degm2 dNcgts_upper_degm2'
        s.write('%s\n'%header)
        print header
        for ibin in range(nbins-1):
            ccc=counts[ibin:].sum()*CORR_BINS[ibin]/SURVEY_AREA
            line='%f %f %f %i %e %e %e %e %f %i %i %i' % \
              (bins[ibin],bins[ibin+1],median_bins[ibin],int(round(counts[ibin])),\
               dn_by_ds[ibin]/(sqDeg2sr*SURVEY_AREA),\
               dn_by_ds_eucl[ibin]/(sqDeg2sr*SURVEY_AREA),\
               dn_by_ds_errs[ibin]/(sqDeg2sr*SURVEY_AREA),\
               dn_by_ds_errs[ibin]/(sqDeg2sr*SURVEY_AREA),\
               CORR_BINS[ibin],\
               round(ccc),round(sqrt(ccc)),round(sqrt(ccc)))
            s.write('%s\n'%line)
            print line
        print '-> There are %i sources in total' % counts.sum()
        s.close()

        print '-> Look in %s' % outputf

        # Now plot a histogram of fluxes to file, with fine binning
        print 'Flux range/uJy = %f -> %f' % (cat[:,BIN_COL].min(),cat[:,BIN_COL].max())
        fig = plt.figure()
        #nHistoBins=100
        binwidth=50.0
        #n,b,p=plt.hist(cat[:,BIN_COL],nHistoBins,histtype='step',color='black')
        n,b,p=plt.hist(cat[:,BIN_COL], bins=numpy.arange(bins[0],(20.0*SURVEY_NOISE)+binwidth,binwidth),histtype='step',color='black')
        plt.yscale('log')
        plt.xlim(bins[0],20.0*SURVEY_NOISE)
        plt.ylim(0.5,1.0e3)
        plt.xlabel('S/$\mu$Jy')
        plt.ylabel('Number of objects')
        y = numpy.max(n)*gaussian(b,0.0,SURVEY_NOISE,norm=False)
        plt.plot(b,y,'r--',linewidth=1)
        plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2)
        plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
        #plt.text(SURVEY_NOISE,0.16,'1 sigma',rotation=90,color='b',alpha=0.5)
        #plt.text(5.0*SURVEY_NOISE,0.16,'5 sigma',rotation=90,color='b',alpha=0.5)

        #plt.title(r'$\mathrm{Histogram\ of\ IQ:}\ \mu=100,\ \sigma=15$')
        plt.savefig(BOUT_HISTO)
        print '-> Look in %s' % BOUT_HISTO


    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

