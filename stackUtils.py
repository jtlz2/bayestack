"""
Support functions for extractor, stacker and video_recorder
Based on stacker.py

Jonathan Zwart
May 2015

"""

import numpy
import pyfits
from profile_support import profile

#-------------------------------------------------------------------------------

@profile
def readRadio(fradio=None):
    """
    Read radio image into memory (whether map or noise map)
    """

    fits = pyfits.open(fradio)
    imdata=fits[0].data
    imhead=fits[0].header
    fits.close()

    print 'Read map %s' % fradio
    
    return imdata,imhead

#-------------------------------------------------------------------------------

def secateur(incat,BINCOL,cutsDict,numNoiseZone):
    """
    Make cuts on the catalogue here
    Excise stars - 0, 1, -1
    Select AGNn - 0, 1, 2
    Select XXX
    Split into noise zones by noise bin value
    Dictionary format: label:[flag,column]
    """

    for cut,cval in cutsDict.iteritems():
        if cut in ['star','lacy','stern','donley']:
            [ccol,cnum]=cval
            #print cut,cnum,ccol
            if cnum < 0: # i.e. no cut
                incat=incat
            else: # leave 0 or 1
                expr=numpy.where(incat[:,ccol].astype('int')==cnum)
                incat=incat[expr]
        elif 'noise%i'%numNoiseZone in cut:
            [ccol,clow,chigh]=cval # clow < noise < chigh
            noises=numpy.power(incat[:,ccol],-0.5) # w -> rms
            expr=numpy.logical_and(noises>clow,noises<chigh)
            incat=incat[expr]

    return incat

#-------------------------------------------------------------------------------

