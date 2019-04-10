#!/usr/bin/env python

"""
Resample a radio map and write to file

***BEWARE: No header info or WCS is written, so these files should not
(yet) be used for e.g. stacking
i.e. Better to do this on-the-fly inside extractor.py/stacker.py
"""

import sys
import importlib
import numpy,pyfits
import stacker
from profile_support import profile
from image_registration.fft_tools import upsample_image
from image_registration.fft_tools.downsample import downsample
from utils import block_mean

__name_cached=__name__
param_file=sys.argv[-1]

#setf='%s' % param_file.split('.')[-2]
# Import the settings variables
if __name__=='__main__':
    setf='bayestack_settings'
else:
    setf='%s'%param_file

set_module=importlib.import_module(setf)
globals().update(set_module.__dict__)
__name__=__name_cached

#-------------------------------------------------------------------------------

@profile
def main():

    """
    """

    CLOBBER=True
    # Read the radio data
    imdata,imhead=stacker.readRadio(fradio=inputMap)
    # Read the noise map
    noisedata,noisehead=stacker.readRadio(fradio=inputNoiseMap)

    # Resample
    use_numpy_fft=True
    print '--> Upsampling radio map x %i' % upsamplingFactor
    print      '%s -> %s' % (inputMap,outputMap)

    # Handle NaN blanking (if present, otherwise return the array)
    # Adding a large offset (or NaNs) causes dynamic range issues (striping), so set NaN -> 0 instead
    #im[numpy.isnan(im)]=0.0#-999.0#-imhead['DATAMIN']
    # Blanked pixels are now zeros, which suits the upsampler
    im=numpy.nan_to_num(imdata[0,0,:,:])
    zmask = (im==0)

    upsampledMap=upsample_image(im,\
                                    upsample_factor=upsamplingFactor,\
                                    use_numpy_fft=use_numpy_fft)

    # Reapply the original mask:
    #http://stackoverflow.com/questions/32846846/quick-way-to-upsample-numpy-array-by-nearest-neighbor-tiling
    #upsampledMap[upsampledMap<-998.0]=0.0#numpy.NaN
    upsampledMask=zmask.repeat(upsamplingFactor,axis=0).repeat(upsamplingFactor,axis=1)
    upsampledMap=numpy.ma.masked_array(upsampledMap,mask=upsampledMask).filled(fill_value=0.0)

    # Fix header keywords for new pixel scale
    for p in [1,2]:
        imhead['NAXIS%i'%p] *=    int(upsamplingFactor)
        imhead['CDELT%i'%p] *=    1.0/abs(upsamplingFactor)
        imhead['CRPIX%i'%p]  =    float(int(imhead['NAXIS%i'%p])/2 +1)
        noisehead['NAXIS%i'%p] *= int(upsamplingFactor)
        noisehead['CRPIX%i'%p]  = float(int(noisehead['NAXIS%i'%p])/2 +1)
        noisehead['CDELT%i'%p] *= 1.0/abs(upsamplingFactor)

    pyfits.writeto(outputMap,upsampledMap,header=imhead,output_verify='fix',clobber=CLOBBER)
    print '--> Finished upsampling radio map'
    #sys.exit(0)

    # Calculate a residual image and generate some stats
    if False:
        print im.shape
        x=block_mean(upsampledMap,upsamplingFactor)
        #numpy.roll(x,20,axis=0)
        #numpy.roll(x,20,axis=1)
        print x.shape
        print im.size
        print x.size
        print numpy.count_nonzero(~numpy.isnan(im))
        print numpy.count_nonzero(~numpy.isnan(x))
        print numpy.count_nonzero(numpy.isnan(im))
        print numpy.count_nonzero(numpy.isnan(x))

        print numpy.nanstd(im)
        print numpy.nanstd(x)
        print numpy.nanmean(im)
        print numpy.nanmean(x)
        r=im-x
        #pyfits.writeto('r.fits',r,clobber=CLOBBER)
        #numpy.savetxt('r.txt',r.flatten())
        print numpy.nanstd(r)
        print numpy.nanmean(r)
        #print r[~mask].shape
        print r.size
        #print r[~mask].size
        #print numpy.nanstd(r[~mask])
        #print numpy.nanmean(r[~mask])
        y=downsample(upsampledMap,upsamplingFactor)
        q=im-y
        pyfits.writeto('q.fits',q,clobber=CLOBBER)
        #sys.exit(0)

    # Stop here if in P mode [use I mode to handle noise map also]
    if doPoln:
        print 'Skipping upsampling of noise map'
        return 0

    # Upsample the noise map
    print '--> Upsampling noise map x %i' % upsamplingFactor
    print '    %s -> %s' % (inputNoiseMap,outputNoiseMap)

    # Handle NaN blanking
    n=numpy.nan_to_num(noisedata[0,0,:,:])
    #n[numpy.isnan(n)]=-999.0#-noisehead['DATAMIN']
    nmask = (n==0)
    upsampledNoiseMap=upsample_image(n,\
                                    upsample_factor=upsamplingFactor,\
                                    use_numpy_fft=use_numpy_fft)

    # Reapply the original mask:
    upsampledNoiseMask=nmask.repeat(upsamplingFactor,axis=0).repeat(upsamplingFactor,axis=1)
    upsampledNoiseMap=numpy.ma.masked_array(upsampledNoiseMap,mask=upsampledNoiseMask).filled(fill_value=0.0)
    #upsampledNoiseMap[upsampledNoiseMap<-998.0]=numpy.NaN

    pyfits.writeto(outputNoiseMap,upsampledNoiseMap,output_verify='fix',clobber=CLOBBER)
    print '--> Finished upsampling noise map'

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)


    
