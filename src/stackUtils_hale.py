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

def cmask(index,radius,array):
    a,b = index
    nx,ny = array.shape
    y,x = np.ogrid[-a:nx-a,-b:ny-b]
    mask = x*x + y*y <= radius*radius

    return array[mask]

#-----------------------------------------------------------

def tophat(height,size,radius):
    """
    """
    x, y = numpy.mgrid[-size//2 + 1:size//2 + 1, -size//2 + 1:size//2 + 1]
    mask = x*x + y*y <= radius*radius
    t = numpy.zeros((size,size))
    t[mask] = height

    return t,mask

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

def calculateNoiseZones(wvlaf,noiseRanges,SURVEY_NOISE,noisezonesf,maskf=None):
    """
    Given a *sensitivity* map
    and a list of noise ranges, calculate the number of pixels in each range
    """

    uu=pyfits.open(wvlaf)
    try:
        W_VLA=uu[0].data[0,0,:,:]
    except IndexError:
        W_VLA=uu[0].data[:,:]
    D1 = uu[0].header['CDELT1']
    D2 = uu[0].header['CDELT2']
    uu.close()
    SAperPix = abs(D1*D2) # Area per pixel in *sq. deg.*
    #rmsmap=SURVEY_NOISE*numpy.power(W_VLA[W_VLA!=0.0],-0.5)
    #rmsmap=SURVEY_NOISE*numpy.power(W_VLA,-0.5)
    rmsmap=W_VLA

    # Open up the mask
    if maskf is not None:
        mm=pyfits.open(maskf)
        try:
            M_VLA=mm[0].data[0,0,:,:]
        except IndexError:
            M_VLA=mm[0].data[:,:]
        mm.close()
        mask=M_VLA

        pixels=rmsmap[mask==0.0].flatten()
    #print pixels.min(),pixels.max()
    else:
        pixels=rmsmap.flatten()

    nz=open(noisezonesf,'w')
    hdr='# zone noise_min_uJy noise_max_uJy area_sq_deg'
    print hdr
    nz.write('%s\n'%hdr)
    noiseAreas=numpy.zeros(len(noiseRanges))
    for i,[noiseMin,noiseMax] in enumerate(noiseRanges):
        expr=numpy.logical_and(pixels>noiseMin,pixels<noiseMax)
        noiseAreas[i]=SAperPix*pixels[expr].size
        line='%i %6.4f %6.4f %e' % (i,noiseMin,noiseMax,noiseAreas[i])
        print line
        nz.write('%s\n'%line)
        #print noiseMin,noiseMax,noiseAreas
    nz.close()
    print 'Look in %s' % noisezonesf

    del W_VLA,rmsmap

    return noiseAreas

#-------------------------------------------------------------------------------

def secateur(incat,cutsDict,numNoiseZone):
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
            #noises=SURVEY_NOISE*numpy.power(incat[:,ccol],-0.5) # w->rms
            noises=incat[:,ccol]
            print noises.min(),noises.max()
            expr=numpy.logical_and(noises>clow,noises<chigh)
            incat=incat[expr]

    return incat

#-------------------------------------------------------------------------------

