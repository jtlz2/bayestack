#!/usr/bin/env python

"""
This is picker.py
Jonathan Zwart
November 2015

Needs to be integrated into main set of scripts

Usage:

NB randomPosns=False/True
Set catf below
Set upsamplingFactor below
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_150520.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_servsmasked_150813.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_random_posns_150813.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_servsmasked_w_poln_150818.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_servsmasked_shift10p10p_w_poln_151001.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_servsmasked_random_posns_w_poln_151001.txt

./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_servsmasked_w_poln_151116.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_servsmasked_shift800p800p_w_poln_151116.txt
./picker.py > /Users/jtlz2/elaisn1/servs/vla_overlap_servsmasked_random_posns_w_poln_151116.txt

Then cross match with STILTS:
match-vla.stil
match-with-vla-detns.stil

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

    upsamplingFactor=8
    # Read radio I map
    radiof='/Users/jtlz2/elaisn1/maps/EN1_WITHBEAM_vla_x%i.fits'%abs(upsamplingFactor)
    data,hdr=readRadio(fradio=radiof)
    vla_map=data[:,:]
    #print 'Read %s' % radiof
    WCS=astWCS.WCS(hdr,mode='pyfits')
    xmax=hdr['NAXIS1']
    ymax=hdr['NAXIS2']

    # Read radio P map
    polf='/Users/jtlz2/elaisn1/maps/EN1.Pave.mosaic_vla_x%i.fits'%abs(upsamplingFactor)
    pol,phdr=readRadio(fradio=polf)
    vla_polmap=pol[:,:]
    #print 'Read %s' % radiof
    WCSpol=astWCS.WCS(phdr,mode='pyfits')
    xmaxp=phdr['NAXIS1']
    ymaxp=phdr['NAXIS2']

    # Read radio sensitivity map
    weightf='/Users/jtlz2/elaisn1/maps/EN1.I.mosaic.sensitivity_vla_x%i.fits'% abs(upsamplingFactor)
    weight,whdr=readRadio(fradio=weightf)
    weight_map=weight[:,:]

    # Read catalogue
    catf='/Users/jtlz2/elaisn1/servs/servs-en1-full-data-fusion-sextractor-cutdown-masked-150522.txt'
    cat=numpy.genfromtxt(catf)
    ngals=cat.shape[0]

    randomPosns=False
    posnShiftX=0.0 # pixels, or 0.0
    posnShiftY=0.0 # pixels, or 0.0
    if randomPosns:
        ngals=15779#15226
        numpy.random.seed(seed=1234)
    ncumul=0
    print '# n id RA Dec x_radio y_radio S_4p8 SoverSqrtW_4p8 W_4p8 P_4p8 PoverSqrtW_4p8'
    for igal in xrange(ngals):
        idd=cat[igal,0]; ra=cat[igal,1]; dec=cat[igal,2]
        if randomPosns:
            notOK=True
            while(notOK):
                x=numpy.random.randint(0,hdr['NAXIS1'])
                y=numpy.random.randint(0,hdr['NAXIS2'])
                if vla_map[y,x]!=0: # This condition needs to be revisited after upsampling
                    notOK=False
                    ncumul+=1
                    print '%i %i %f %f %f %f %e %e %e %e %e'\
                           %(ncumul,idd,0.0,0.0,x,y,vla_map[y,x],\
                             vla_map[y,x]/numpy.power(weight_map[y,x],0.5),\
                             weight_map[y,x],vla_polmap[y,x],\
                             vla_polmap[y,x]/numpy.power(weight_map[y,x],0.5))

        else:
            x,y=astWCS.WCS.wcs2pix(WCS,ra,dec)
            x+=posnShiftX; y+=posnShiftY
            if x<0 or x>xmax or y<0 or y>ymax or vla_map[y,x]==0:
                continue
            else:
                ncumul+=1
                print '%i %i %f %f %f %f %e %e %e %e %e'\
                  %(ncumul,idd,ra,dec,x,y,vla_map[y,x],\
                    vla_map[y,x]/numpy.power(weight_map[y,x],0.5),\
                    weight_map[y,x],vla_polmap[y,x],\
                    vla_polmap[y,x]/numpy.power(weight_map[y,x],0.5))

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
