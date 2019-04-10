#!/usr/bin/env python

"""
Based on Tessa Vernstrom's code
22 August 2014
"""

import os,sys
import numpy
from numpy import exp,log10
from scipy import fftpack

#from scipy import *
#import pylab as pyl
#from scipy import stats

#-------------------------------------------------------------------------------

def setupDArray(ss):
    """
    #SETUP D ARRAY
    """
    dnn = 2**18 #SIZE OF D ARRAY SOMETHING LARGE FOR GOOD RESOLUTION

    dmin = 0. #0 IF NO NOISE
    dmax = ss.max() #CAN CHANGE THIS BUT LIKELY WANT AT LEAST THE MAXIMUM VALUE OF THE FLUX DENSITY SOURCE COUNT YOU INPUT
    dx  =(dmax-dmin)/(dnn-1.)
    dd = numpy.arange(dmin,dmax+dx/2.,dx)
    dd2 = (dd[1:]+dd[:-1])/2.
    #dnn=len(dd)
    return dd,dx



#-------------------------------------------------------------------------------

def setupBeam(b_fwhm,pxsz_as):
    """
    #SETUP BEAM ARRAY (ASSUME GAUSSIAN)
    """


    b_fwhmpx = b_fwhm/pxsz_as#BEAM FWHM IN PIXELS
    bx = numpy.arange(b_fwhmpx*3.) #CREATE X ARRAY TO CALC GAUSSIAN, PROBABLY 2 OR 3 TIMES THE SIZE
    bx = bx-(len(bx)-1)/2
    beam = exp((-1)*((bx.reshape(1,len(bx))**2)+(bx.reshape(len(bx),1)**2))/(2.*(b_fwhmpx/2.3548)**2))
    beam1 = beam.reshape(beam.size)
    beam1=beam1[beam1>=1.e-06] #SET SOME LOWER LIMIT
    nb=len(beam1)
    return beam1

#-------------------------------------------------------------------------------


def generateRx(beam1,pxsz_as,ss,dnds,dd):
    """
    #DO RX INTEGRAL
    """

    dnn=len(dd)
    logs=log10(ss)

    #PIXEL SIZE IN STERADIANS
    pxsz_sr = pxsz_as/4.25452e10
    #SCALE COUNTS TO PER PIXEL INSTEAD OF PER SR
    logdnds_px = log10(dnds)+log10(pxsz_sr)

    rx = numpy.zeros(dnn)
    for i in range(dnn):
        flux = dd[i]/beam1
        flux1=flux[(flux>=ss.min())&(flux<ss.max())]
        beam2=beam1[(flux>=ss.min())&(flux<ss.max())]
        logflx=log10(flux1)
        newy=numpy.interp(logflx,logs,logdnds_px)
        rx[i]=sum(10.**newy/beam2)

    return rx


#-------------------------------------------------------------------------------

def fftRx(rx,dx,sigma=None,dnn=None):
    """
    """
    #COMPUTE FFTS
    frx0=rx*dx
    frx1=fftpack.fft(frx0)
    frx2=frx1-frx1[0]

    ##IF YOU WANT TO ADD NOISE
    #sigma=
    #omega=fftpack.fftfreq(dnn,dx)*numpy.pi
    #frx2=frx2-(sigma**2*omega**2)/2.

    frx3=exp(frx2)

    #INVERSE IT TO GET P(D), THIS SHOULD SUM TO 1
    pofd=numpy.real(fftpack.ifft(frx3))
    return pofd

#-------------------------------------------------------------------------------

def pofd2confNoise(dd,pofd):
    """
    #COMPUTE CONFUSION NOISE    
    """
    pofd_csum=numpy.cumsum(pofd)
    ll=dd[numpy.argmin(abs(pofd_csum-.16))]
    uu=dd[numpy.argmin(abs(pofd_csum-.84))]
    conf=(uu-ll)/2.
    conf=conf/1.e-06

    return conf

#-------------------------------------------------------------------------------

def dnds2conf(ss,dnds,b_fwhm,pxsz_as):
    """
    Convert a source count to an effective confusion noise given a beam
    """

    dd,dx=setupDArray(ss)
    beam1=setupBeam(b_fwhm,pxsz_as)
    rx=generateRx(beam1,pxsz_as,ss,dnds,dd)
    pofd=fftRx(rx,dx)
    conf=pofd2confNoise(dd,pofd)
    print 'confusion noise=',conf.round(3),' microJy/beam'

    return conf

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # These are the data
    #LOG(S), DNDS, LOG(DNDS)
    print '***NB Example:***'
    ss = 10.0**numpy.array([-7.32,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9])
    dnds = 10.0**numpy.array([16.1737,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314]) #log10(dnds)

    b_fwhm = 8. #BEAM FWHM IN ARCSECONDS
    pxsz_as = 1.5667 #PIXEL SIZE IN ARCSECONDS
    conf=dnds2conf(ss,dnds,b_fwhm,pxsz_as)
    print conf

#-------------------------------------------------------------------------------

    sys.exit(0)

