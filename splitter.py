#!/usr/bin/env python

"""
This is splitter.py
Jonathan Zwart
June 2015

Make catalogue cuts and
Bin data set for bayestack.py
Based on binner.py

Usage:

./splitter.py SETTINGS_FILE.py

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
import stackUtils
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
        cat[:,BIN_COL] *= 1000.0 # mJy -> uJy in SURVEY_AREA sq. deg.
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

    cutsDict={'star':[30,0],'lacy':[34,-1],'stern':[38,-1],'donley':[40,-1],\
              'noise1':[47,1.0,2.5],'noise2':[47,2.5,5.0],\
              'noise3':[47,5.0,8.0],'noise4':[47,8.0,12.0]}

    numNoiseZones=len([k for k in cutsDict.keys() if 'noise' in k])

    # Need to calculate survey areas

    # Set up the plot
    fig = plt.figure()
    colors=['y','g','b','k']

    idl_s=False
    for n in range(numNoiseZones):
        f='.'.join(['_'.join([BOUT_CAT.split('.')[-2],'a%i'%n]),'txt'])
        print f,n
        print cat.shape
        ccat=stackUtils.secateur(cat,BIN_COL,cutsDict,n)
        print ccat.shape
        countUtils.writeCountsFile(f,bins,ccat,SURVEY_AREA,\
                               idl_style=idl_s,verbose=False,corrs=CORR_BINS)
        binwidth=0.1*SURVEY_NOISE
        #sys.exit(0)
        n,b,p=plt.hist(ccat[:,BIN_COL], bins=numpy.arange(bins[0],(20.0*SURVEY_NOISE)+binwidth,binwidth),histtype='step',color=colors[n])

    if True:
        # Now plot a histogram of fluxes to file, with fine binning
        print 'Flux range/uJy = %f -> %f' % (ccat[:,BIN_COL].min(),ccat[:,BIN_COL].max())
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
        #plt.title('')
        plt.savefig(BOUT_HISTO)
        print '-> Look in %s' % BOUT_HISTO

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

