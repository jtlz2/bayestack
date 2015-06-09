"""
Support functions for extractor, stacker and video_recorder
Based on stacker.py

Jonathan Zwart
May 2015

"""

import os,sys
import importlib,glob
import numpy
import pyfits
from profile_support import profile

if 'chains' in sys.argv[-1]:
    potential_settings=glob.glob(os.path.join(sys.argv[-1],'*settings*py'))
    assert(len(potential_settings)==1), '***More than one potential settings file!'
    settingsf='.'.join([sys.argv[-1],potential_settings[0].split('/')[-1].split('.')[-2]])
else:
    settingsf=sys.argv[-1].split('.')[-2]

print '%s is using %s' % (__name__,settingsf)
try:
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)
except:
    print '***Warning: Settings not loaded'


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
            noises=SURVEY_NOISE*numpy.power(incat[:,ccol],-0.5) # w->rms
            print noises.min(),noises.max()
            expr=numpy.logical_and(noises>clow,noises<chigh)
            incat=incat[expr]

    return incat

#-------------------------------------------------------------------------------

