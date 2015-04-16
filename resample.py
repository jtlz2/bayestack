#!/usr/bin/env python

"""
Resample a radio map and write to file

***BEWARE: No header info or WCS is written, so these files should not
(yet) be used for e.g. stacking
i.e. Better to do this on-the-fly inside extractor.py/stacker.py
"""

import sys
import importlib
import pyfits
import stacker
from profile_support import profile
from image_registration.fft_tools import upsample_image


__name_cached=__name__
param_file=sys.argv[-1]

#setf='%s' % param_file.split('.')[-2]
# Import the settings variables
if __name__=='__main__':
    setf='settings'
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

    # Read the radio data
    imdata,imhead=stacker.readRadio(fradio=inputMap)
    # Read the noise map
    noisedata,noisehead=stacker.readRadio(fradio=inputNoiseMap)

    # Resample
    use_numpy_fft=True
    print '--> Upsampling radio map x %i' % upsamplingFactor
    print      '%s -> %s' % (inputMap,outputMap)
    upsampledMap=upsample_image(imdata[0,0,:,:],\
                                    upsample_factor=upsamplingFactor,\
                                    use_numpy_fft=use_numpy_fft)
    pyfits.writeto(outputMap,upsampledMap)
    print '--> Finished upsampling radio map'

    print '--> Upsampling noise map x %i' % upsamplingFactor
    print '    %s -> %s' % (inputNoiseMap,outputNoiseMap)
    upsampledNoiseMap=upsample_image(noisedata[0,0,:,:],\
                                    upsample_factor=upsamplingFactor,\
                                    use_numpy_fft=use_numpy_fft)
    pyfits.writeto(outputNoiseMap,upsampledNoiseMap)
    print '--> Finished upsampling noise map'

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)


    
