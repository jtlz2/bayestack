#!/usr/bin/env python

"""
This is picker.py
Jonathan Zwart
May 2015

Usage:

NB randomPosns=False/True
Set catf below
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap.txt
picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_150520.txt

Then cross match with STILTS

"""

import os,sys
from astLib import astWCS
import numpy
from profile_support import profile
from stackUtils import readRadio

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file

#-------------------------------------------------------------------------------


@profile
def main():

    """
    """

    # Read radio map
    radiof='/Users/jtlz2/Dropbox/elaisn1/maps/EN1_WITHBEAM_vla.FITS'
    data,hdr=readRadio(fradio=radiof)
    vla_map=data[0,0,:,:]
    #print 'Read %s' % radiof
    WCS=astWCS.WCS(hdr,mode='pyfits')
    xmax=hdr['NAXIS1']
    ymax=hdr['NAXIS2']

    # Read radio sensitivity map
    weightf='/Users/jtlz2/Dropbox/elaisn1/maps/EN1.I.mosaic.sensitivity_vla.fits'
    weight,whdr=readRadio(fradio=weightf)
    weight_map=weight[0,0,:,:]

    # Read catalogue
    catf='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-150522.txt'
    cat=numpy.genfromtxt(catf)
    ngals=cat.shape[0]

    randomPosns=False
    if randomPosns:
        ngals=15226
        numpy.random.seed(seed=1234)
    ncumul=0
    print '# n id RA Dec x_radio y_radio S_4p8 SoverW_4p8 W_4p8'
    for igal in xrange(ngals):
        idd=cat[igal,0]; ra=cat[igal,1]; dec=cat[igal,2]
        if randomPosns:
            notOK=True
            while(notOK):
                x=numpy.random.randint(0,hdr['NAXIS1'])
                y=numpy.random.randint(0,hdr['NAXIS2'])
                if vla_map[y,x]!=0:
                    notOK=False
                    ncumul+=1
                    print '%i %i %f %f %f %f %e %e %e'\
                           %(ncumul,idd,ra,dec,x,y,vla_map[y,x],\
                             vla_map[y,x]/weight_map[y,x],weight_map[y,x])
        else:
            x,y=astWCS.WCS.wcs2pix(WCS,ra,dec)
            if x<0 or x>xmax or y<0 or y>ymax or vla_map[y,x]==0:
                continue
            else:
                ncumul+=1
                print '%i %i %f %f %f %f %e %e %e'\
                  %(ncumul,idd,ra,dec,x,y,vla_map[y,x],vla_map[y,x]/weight_map[y,x],\
                    weight_map[y,x])

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
