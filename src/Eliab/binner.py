#!/usr/bin/env python

"""
This is binner.py
Jonathan Zwart
June 2015

Bin data set for bayestack.py

Usage:

./binner.py SETTINGS_FILE.py

BIN_CAT_FORMs:
0 COSMOS-VLA data [Ketron]
1 VIDEO catalogue from stacker.py [JZ]
2 VVDF Bondi 2003 catalogue from web [web]
3 LR-matched Bondi 2003 catalogue [Kim]
4 VIDEO catalogue from stacker.py [JZ]
5 10C_LH catalogue [Imogen]
6 SDSS catalogue [Eliab]
7 10C_LH tier 2 [Imogen]
8 ELAIS-N1 JVLA [JZ]
9 ELAIS-N1 GMRT [JZ]

*** IF YOUR CATALOGUE IS IN MUJY YOU NEED
    TAKE NO FURTHER ACTION IN THIS FILE ***

"""

import os,sys
import importlib
import numpy
from profile_support import profile
from utils import *
import countUtils
import matplotlib.pyplot as plt

param_file=sys.argv[-1]
setf='%s' % param_file.split('.')[-2]
print '%s is using %s' % (__name__,setf)

#-------------------------------------------------------------------------------

@profile
def main():
    """
    """

    global CORR_BINS

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    if CORR_BINS is None:
        CORR_BINS=numpy.ones(len(bins)-1)

    print 'Reading from %s' % BIN_CAT
    cat=numpy.genfromtxt(BIN_CAT)

    # Convert unit if required
    if BIN_CAT_FORM in [0,2,3]:
        cat[:,BIN_COL] *= 1e6#1000.0 # mJy -> uJy in SURVEY_AREA sq. deg.
        # 5 sigma:
        #cat=cat[numpy.where((cat[:,BIN_COL]/cat[:,BIN_COL+1])>0.0)]
    elif BIN_CAT_FORM in [6,8]:
        cat[:,BIN_COL] *= Jy2muJy
    else:
        pass

    # Check the corrections
    assert(len(CORR_BINS) == len(bins)-1),\
      '**Binning corrections mismatch %s' % (BIN_CAT,bins,CORR_BINS)

    # Correct the fluxes for the resolution bias
    if CORR_RESOLUTION is not None:
        print '--> Corrected fluxes for resolution bias (x %f)' % CORR_RESOLUTION
        cat[:,BIN_COL] *= CORR_RESOLUTION

    # Optionally threshold the catalogue
    if BIN_CAT_CLIP is not None:
        cat=cat[numpy.where(cat[:,BIN_COL_CLIP]>BIN_CAT_CLIP)]
        #cat=cat[numpy.where(numpy.abs(cat[:,BIN_COL_CLIP])>BIN_CAT_CLIP)]

    #print 'S/uJy: %f -> %f' % (numpy.min(cat[:,BIN_COL]),numpy.max(cat[:,BIN_COL]))

    idl_s=False
    countUtils.writeCountsFile(BOUT_CAT,bins,cat[:,BIN_COL],SURVEY_AREA,\
                               idl_style=idl_s,verbose=True,corrs=CORR_BINS,redshift=cat[:,2])

    if True:
        # Now plot a histogram of fluxes to file, with fine binning
        zbins =[0.2, 0.45, 0.7, 1.0, 1.3, 1.6, 1.85, 2.15, 2.35, 2.55, 2.85, 3.15, 3.5]
        z_min = zbins[int(num) -1]
        z_max = zbins[int(num)]
        print 'Flux range/uJy = %f -> %f' % (cat[:,BIN_COL].min(),cat[:,BIN_COL].max())
        fig = plt.figure()
        binwidth=0.2*SURVEY_NOISE
        b=plt.hist(cat[:,BIN_COL], bins=numpy.arange(bins[0],(20.0*SURVEY_NOISE)+binwidth,binwidth),histtype='step',color='black')[1]
        n = len(cat)
        plt.yscale('log')
        plt.xlim(bins[0],20.0*SURVEY_NOISE)
        plt.ylim(0.5,3.0e3)
        plt.xlabel('S/$\mu$Jy',fontsize=18)
        plt.ylabel('Number of objects',fontsize=18)
        y = 0.06*numpy.max(n)*gaussian(b,80.0,SURVEY_NOISE,norm=False)
        plt.plot(b,y,'r--',linewidth=1)
        plt.axvline(1.0*SURVEY_NOISE,color='g',alpha=0.2)
        plt.axvline(5.0*SURVEY_NOISE,color='g',alpha=0.2)
        plt.text(SURVEY_NOISE,1e3,'$ \sigma$',rotation=90,color='b',alpha=0.5,fontsize=18)
        plt.text(5.0*SURVEY_NOISE,1e3,'$5 \sigma$',rotation=90,color='b',alpha=0.5,fontsize=18)
        plt.text(1.5e3,1e3,'%.2f < z < %.2f\n N = %d'%(z_min,z_max,len(cat)),alpha=0.5,fontsize=18 )
        #plt.title('')
        plt.savefig(BOUT_HISTO)
        print '-> Look in %s' % BOUT_HISTO

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

