"""
Support functions for extractor, stacker and video_recorder
Based on stacker.py

Jonathan Zwart
May 2015

"""


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


