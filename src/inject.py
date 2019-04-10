#!/usr/bin/env python

"""

readcol,'rqq_sfr10.cat',clust,sftype,agntype,ra,dec,dist1,z,z_mod,itot_151,itot_610,itot_1400,format='x,f,f,f,d,d,d,d,f,f,f,f'


;stop

resolution=5.   ; arcsec
fluxlimit=1.0d-5 ; Jansky/beam

npixpix=38000

npix=npixpix/resolution
print,'npix=',npix
arr=fltarr(npix,npix)

inp=readfits('Image3.fits',h1)

adxy,h1,ra,dec,x1,y1

flagra=fltarr(n_elements(ra))
flagdec=fltarr(n_elements(ra))
w1=where(ra ge 0.)
w2=where(ra lt 0.)
w3=where(dec ge 0.)
w4=where(dec lt 0.)
flagra=npixpix/2./resolution
;flagra(w2)=1900.
flagdec=npixpix/2./resolution
;flagdec(w4)=1900.


x1=ra*(-36000./resolution)+flagra
y1=dec*(36000./resolution) +flagdec


rmsflux=fluxlimit/ sqrt(!pi*(0.5*resolution)^2.)

;print,'npix=',floor(npix)

;print,'ran=',randomn(seed,npix,npix)*rmsflux
arr=randomn(seed,npix,npix)*rmsflux



sxaddpar,h1,'CRPIX1',(npixpix+1000.)/2./resolution 
sxaddpar,h1,'CRPIX2',(npixpix+1000.)/2./resolution
sxaddpar,h1,'CDELT1',(-1./3600.)/(10./resolution)
sxaddpar,h1,'CDELT2',(1./3600.)/(10./resolution)
sxaddpar,h1,'NAXIS',2.

sxaddpar,h1,'CRVAL1',180.
sxaddpar,h1,'CRVAL2',0.
sxaddpar,h1,'BUNIT','JY/BEAM'


;psf=psf_gaussian(npixel=51,fwhm=[resolution,resolution])


;for i=0,3599 do begin
;

imconv=convolve(arr,psf)
writefits,'kakcon.fits',imconv,h1
for j=0,n_elements(ra)-1 do begin
;print,ra(j),dec(j),x1(j),y1(j)
psfsource=psf_gaussian(npixel=101,fwhm=[resolution,resolution])*10.^(itot_1400(j))
xmin=round(x1(j)-50)
xmax=round(x1(j)+50)
ymin=round(y1(j)-50)
ymax=round(y1(j)+50)

;print,ra(j)+180.,dec(j)
cut=imconv[xmin:xmax,ymin:ymax]
imconv[xmin:xmax,ymin:ymax]=cut+psfsource
endfor
;endfor


;writefits,'kak.fits',arr,h1
writefits,'Map_5arcsec_10uJy.fits',imconv,h1

end


Resample a radio map and write to file
"""

import os,sys
import importlib
import pyfits
import stacker
from profile_support import profile
import numpy,scipy
from numpy.random import RandomState
from math import pi,log,erf,sqrt
#if os.getenv('PBS_O_HOST')==None or os.getenv('PBS_O_HOST').split('.')[0]!='baltasar':
#    from photutils.psf import GaussianPSF
#    from photutils import psf_photometry
from image_registration.fft_tools import upsample_image
from utils import Jy2muJy
from utils import matchit

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

def gauss_kern(size,sizey=None,sigma=None,norm=False):
    """
    Returns a normalized 2D gauss kernel array for convolutions
    Adapted from http://astrolitterbox.blogspot.com/2012/04/creating-discrete-gaussian-kernel-with.html
    """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)

    if sigma is None:
        sigma=size

    x, y = numpy.mgrid[-size:size+1, -sizey:sizey+1]
    g = numpy.zeros((2*size+1,2*sizey+1))
    g += numpy.exp(-((x+size)**2/(2.0*float(sigma)**2)+(y+sizey)**2/(2.0*float(sigma)**2)))
    if norm: g = 1.0/g.sum()
    return g


#-------------------------------------------------------------------------------

def gaussian2(size, sigma, rpix=None,xshift=0.0,yshift=0.0):
    """Returns a normalized circularly symmetric 2D gauss kernel array
    f(x,y) = A.e^{-(x^2/2*sigma^2 + y^2/2*sigma^2)} where
    A = 1/(2*pi*sigma^2)
    as defined by Wolfram Mathworld
    http://mathworld.wolfram.com/GaussianFunction.html
    Adapted from
    https://github.com/mubeta06/python/blob/master/signal_processing/sp/gauss.py
    Also returns a mask of radius rpix for the gaussian kernel
    """

    A = 1.0#/(2.0*numpy.pi*sigma**2)
    x, y = numpy.mgrid[-size//2 + 1:size//2 + 1, -size//2 + 1:size//2 + 1]
    g = A*numpy.exp(-(((x+xshift)**2/(2.0*sigma**2))+((y+yshift)**2/(2.0*sigma**2))))

    if rpix is not None:
        assert(rpix<=size//2+1), 'rpix too large, %g > %g' % (rpix,size//2+1)
        mask = (x**2 + y**2 <= rpix**2)
        return g,mask
    else:
        return g

#-------------------------------------------------------------------------------

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
def main():

    """
    """

    global upsamplingFactor

    # Read the radio data
    imdata,imhead=stacker.readRadio(fradio=injectInputMap)
    # Read the noise map
    #noisedata,noisehead=stacker.readRadio(fradio=injectInputNoiseMap)

    #upsamplingFactor=8

    if upsamplingFactor is None: upsamplingFactor=1

    if injectFakeNoise:
        print 'Creating new map with fake noise x%i' % upsamplingFactor
        numpy.random.seed(seed=SEED_SIM)
        # Strictly I think -2 and -1 should be reversed
        u=numpy.random.normal(0.0,injectionNoise/Jy2muJy,size=(upsamplingFactor*imdata.shape[-2],upsamplingFactor*imdata.shape[-1]))
        #imdata=numpy.random.normal(0.0,NOISE_SIM/Jy2muJy,size=imdata.shape)
        #from astropy.convolution import Gaussian2DKernel,convolve,convolve_fft
        fwhm=radioSynthBeamFWHM*upsamplingFactor # pixels
        stddev=fwhm/2.35482
        #gauss = Gaussian2DKernel(stddev=stddev)
        #print 'Convolving noise with PSF (stddev=%f)' % stddev
        #u = convolve_fft(u,gauss,allow_huge=True)#,boundary='extend')
        u=scipy.ndimage.filters.gaussian_filter(u,stddev)
    else:
        print 'Creating new map without fake noise'
        # Strictly I think -2 and -1 should be reversed
        u=numpy.zeros((upsamplingFactor*imdata.shape[-2],upsamplingFactor*imdata.shape[-1]))

    #u=upsample_image(imdata[0,0,:,:],upsample_factor=upsamplingFactor,\
    #                 use_numpy_fft=False)

    #print imhead
    #print noisehead

    # Read the fluxes
    #nsrc=-1
    dumpf=DUMP
    #noisydumpf='%s_noisy.txt' % dumpf.split('.')[0]
    fluxesf=os.path.join(outdir,dumpf)
    S=numpy.loadtxt(fluxesf)[:]
    nsrc=len(S)
    if injectNumberOfSources is not None:
        nsrc=injectNumberOfSources
    #print S
    print 'Injecting %i sources (~ avoiding collisions)' % nsrc

    seed=SEED_SIM
    prng=RandomState(seed)

    # Randomize order of the fluxes from dumpf
    print 'Randomizing ordering of fluxes'
    numpy.random.shuffle(S)

    gaussianDim=int(injectPostageStampRadius*upsamplingFactor)
    gaussianDim2=2*gaussianDim+1

    # Generate fake positions and dump them to file
    xy=[]
    x=numpy.zeros(nsrc); y=numpy.zeros(nsrc)

    outposnf=positionsf
    outposnf=os.path.join(outdir,outposnf)

    if injectRecyclePositions:
        posns=numpy.genfromtxt(outposnf,dtype='int')
        x=posns[:,1]; y=posns[:,2]
    elif injectUseSKADSPosns:
        posnsf='heywood/1sqdeg_0p02uJy_w_posns.txt'
        #posnsf='heywood/1sqdeg_0p02uJy_w_posns_sorted_by_flux.txt' # high2low
        skads=numpy.genfromtxt(posnsf)
        #print u.shape
        xs=upsamplingFactor*imdata.shape[-2]*(skads[:,0]+0.5) # 1sq.deg.
        ys=upsamplingFactor*imdata.shape[-1]*(skads[:,1]+0.5) # 1sq.deg.
        #H,xedges,yedges=numpy.histogram2d(RAs,Decs,bins=u.shape)
        #print len(H[numpy.where(H>1)])

        x=[]; y=[]
        p=open(outposnf,'w')
        for isrc in range(nsrc):
            if xs[isrc] > upsamplingFactor*gaussianDim and xs[isrc] < upsamplingFactor*(imhead['NAXIS1']-gaussianDim) and ys[isrc] > upsamplingFactor*gaussianDim and ys[isrc] < upsamplingFactor*(imhead['NAXIS2']-gaussianDim):
            #if xs[isrc] > gaussianDim and xs[isrc] < (imhead['NAXIS1']-gaussianDim) and ys[isrc] > gaussianDim and ys[isrc] < (imhead['NAXIS2']-gaussianDim):
                line='%i %i %i'%(isrc,xs[isrc],ys[isrc])
                p.write('%s\n'%line)
                x.append(xs[isrc]); y.append(ys[isrc])
        p.close()
        nsrc=len(x)
        print '%i sources remain after windowing' % nsrc

    else:
        p=open(outposnf,'w')
        isrc=0
        #for isrc in range(nsrc):
        # Avoid the edges!!
        while (isrc<nsrc):
            xdraw=prng.randint(gaussianDim,imhead['NAXIS1']-gaussianDim)
            ydraw=prng.randint(gaussianDim,imhead['NAXIS2']-gaussianDim)
            if True:
            #if (xdraw,ydraw) not in xy: # ~ avoids collisions...
                x[isrc]=xdraw; y[isrc]=ydraw
                xy.append((xdraw,ydraw))
                line='%i %i %i'%(isrc,x[isrc],y[isrc])
                isrc+=1
                p.write('%s\n'%line)
                #print line
        p.close()

        #grid = numpy.indices(u.shape)
        #idx = zip( grid[0].ravel(), grid[1].ravel() )
        #xy=list(numpy.random.shuffle(idx)[:nsrc])
        #print xy
        #xy=[]; i=0; j=0
        #xdraws=prng.randint(gaussianDim,imhead['NAXIS1']-gaussianDim,nsrc)
        #ydraws=prng.randint(gaussianDim,imhead['NAXIS1']-gaussianDim,nsrc)
        #for x,y in zip(xdraws,ydraws):
        #    if x==y:
        #        xdraws[i]=prng.randint(gaussianDim,imhead['NAXIS1']-gaussianDim,nsrc)
        #    i+=1; j+=1
        #p=open(outposnf,'w')
        #for isrc in range(nsrc):
        #    line='%i %i %i'%(isrc,x[isrc],y[isrc])
        #    p.write('%s\n'%line)
        #p.close()
        
        print 'Fake positions written to %s (non-upsampled)' % outposnf

    if not injectUseSKADSPosns:
        print 'Upsampling positions for internal use'
        x*=upsamplingFactor; y*= upsamplingFactor

    posnNoise=injectPositionNoise
    if injectDoPosnNoise:
        xnoise=numpy.random.normal(0.0,posnNoise*upsamplingFactor,size=nsrc)
        ynoise=numpy.random.normal(0.0,posnNoise*upsamplingFactor,size=nsrc)
    else:
        xnoise=numpy.zeros(nsrc); ynoise=numpy.zeros(nsrc)
    #x[0]=y[0]=upsamplingFactor*2500/2.0
    #xnoisy+=x; ynoisy+=y

    fwhm=radioSynthBeamFWHM*upsamplingFactor # pixels
    sig=fwhm/2.35482

    #g=gauss_kern(gaussianDim,gaussianDim,sigma=fwhm/2.35482)
    #rpix=0.5*upsamplingFactor#8.0#gaussianDim
    # require apercorr/(gmsum/gsum) =~ 1.0
    #for rpix in range(1,gaussianDim):
    #    g,mask=gaussian2(gaussianDim2,sig,rpix=rpix)
    #    gmsum=g[mask].sum()
    #    gsum=g.sum()
    #    apercorr=1.0-numpy.exp(-0.5*rpix*rpix/sig/sig)
    #    print rpix,pi*rpix*rpix,fwhm,sig,gmsum,apercorr,gmsum/gsum,apercorr/(gmsum/gsum)-1.0,1.0-apercorr

    rpix=injectMaskRadius*upsamplingFactor
    #print rpix,sig,rpix/sig,sig/rpix

    print 'Postage stamp is %i x %i upix' % (gaussianDim2,gaussianDim2)
    print 'Aperture mask radius / upix %f' % rpix
    print 'Gaussian FWHM / upix = %f' % fwhm
    print 'Gaussian sigma / upix = %f' % sig
    print 'xshift / pix = %f' % injectGaussianShiftX
    print 'yshift / pix = %f' % injectGaussianShiftY

    #for si in range(33):
    #    s=-4.0+0.25*float(si)
    #    g,mask=gaussian2(gaussianDim2,sig,rpix=rpix,xshift=s,yshift=s)
    #    print s,g.sum(),g[mask].sum(),g[mask].max()

    # Initialize the gaussian aperture and mask
    #print rpix,sig,upsamplingFactor
    #sig=1.69*upsamplingFactor
    print gaussianDim,gaussianDim2,sig,rpix
    g,mask=gaussian2(gaussianDim2,sig,rpix=rpix,\
                     xshift=injectGaussianShiftX*upsamplingFactor,\
                     yshift=injectGaussianShiftY*upsamplingFactor)

    #print s,g.sum(),g[mask].sum()

    #for i in range(2*gaussianDim+1):
    #    for j in range(2*gaussianDim+1):
    #        print i,j,g[i,j]

    # Set up a circular mask
    #rpix=3
    #yy,xx=numpy.ogrid[-gaussianDim2:rpix+1,-rpix:rpix+1]
    #mask = yy**2+xx**2 <= rpix**2
    #print mask

    # For photometry
    #if injectAperturePhot:
    #    PSF=GaussianPSF(float(sig),amplitude=1.0)
    fac=pi/(4.0*log(2.0))
    #norm=g.sum()/fac/fwhm/fwhm

    SEx=-99.0*numpy.ones(nsrc); Spix=numpy.zeros(nsrc)
    pSEx=-99.0*numpy.ones(nsrc); pSpix=numpy.zeros(nsrc)
    Ssp=numpy.zeros(nsrc); Smax=numpy.zeros(nsrc)

    # dump photometry to file
    outphotf=injectPhotometryDumpFile
    outphotf=os.path.join(outdir,outphotf)
    ph=open(outphotf,'w')
    hdr='# n x y S SEx dSEx pSEx Spix dSpix pSpix Ssp Smax meanbSpix\n'
    ph.write(hdr)

    if upsamplingFactor==4: fluxCorrFactor=1.00#1.00932315361
    elif upsamplingFactor==8: fluxCorrFactor=1.00#0.954998930166
    else: fluxCorrFactor=1.0

    # Calculate an effective beam radius from the effective beam area
    beamRadiusEff=sqrt(fac/pi)*fwhm


    # Sort the flux array by descending flux
    if injectSortFluxesDescending:
        print 'Sorting flux array (descending)'
        SS=numpy.sort(S)
        S=SS[::-1]
        numpy.savetxt('S.txt',S)

    # Insert the sources
    #nsrc=68000
    fluxScalingFactor=1.0
    for iS in range(nsrc):
        #S[iS]=1.0e6
        #t,mask=tophat(1.0,gaussianDim2,beamRadiusEff) # sqrt(fac)??
        #print t[mask]
        injectedStamp=g*S[iS]*fluxScalingFactor/Jy2muJy # t->g
        #gg,maskmask=gaussian2(gaussianDim2,sig,rpix=2.35482*sig)
        #injectedStamp=numpy.zeros(numpy.shape(gg))
        #injectedStamp[maskmask]=gg[maskmask]*S[iS]*fluxScalingFactor/Jy2muJy
        #print injectedStamp,injectedStamp[mask]
        #print x[iS],y[iS],fluxScalingFactor,gaussianDim,gaussianDim+1
        u[y[iS]-gaussianDim:y[iS]+gaussianDim+1,x[iS]-gaussianDim:x[iS]+gaussianDim+1]+=injectedStamp

    # Count elements that are close to zero
    #numpy.isclose(u,numpy.zeros(numpy.shape(u)),atol=1e-08)
    # Discount the picture frame
    #v=u[upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS1']-gaussianDim),upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS1']-gaussianDim)]
    #nnz=numpy.count_nonzero(v)
    #ntot=v.size
    #nz=ntot-nnz
    #nz2=v[numpy.where(v==0.0)].size
    #print nz,nz2,ntot
    #assert(nz==nz2)
    #focc=float(nnz)/float(ntot)
    #print ntot,nnz,nz,focc

    apercorr=g.sum()/g[mask].sum()
    disccorr=g.sum()/g[mask].size
    #print apercorr,disccorr
    #print g.sum(),g[mask].sum()
    #print g[mask].size,pi*rpix*rpix
    erf2=1.0-numpy.exp(-0.5*rpix*rpix/sig/sig)
    print 'erf2',erf2,g[mask].sum()/g.sum()
    #afwhm=fac*fwhm*fwhm
    #acirc=pi*rpix*rpix
    #amask=g[mask].size
    #inErf=1.0-erf2 # This is the fraction of the gaussian within the aperture


    # Only once the entire imdata array is populated should photometry
    # be attempted
    #gsum=numpy.sum(g)

    if injectNumberOfSourcesExtracted is not None:
        nsrc=injectNumberOfSourcesExtracted


    # Calculate a random (confusion!?) background flux
    backgroundS=0.0
    if injectDoBackSub is not None:
        backgroundS=injectDoBackSub
        #print 'Calculating a background for subtraction (%i sources)'%nsrc
        #samples=numpy.random.choice(u[upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS1']-gaussianDim),upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS2']-gaussianDim)].flatten(),replace=False,size=nsrc)
        #numpy.savetxt('samples.txt',1.0e6*samples) # Big ~ 1 GB!
        #backgroundS=numpy.mean(samples)
    #print 'Assuming a background mean flux = %f uJy' % 1.0e6*backgroundS

        #s=u[upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS1']-gaussianDim),upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS2']-gaussianDim)].shape
        #coords=matchit(u[upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS1']-gaussianDim),upsamplingFactor*gaussianDim:upsamplingFactor*(imhead['NAXIS2']-gaussianDim)].flatten(),xy,s)
        #numpy.savetxt('coords.txt',coords)
        #tellback=[v[i,j] for i,j in coords]
        #numpy.savetxt('tellback.txt',1.0e6*tellback)
        #assert(numpy.allclose(tellback,xy))
        #backgroundS=tellback
    #    offset=1000
    #    samples=numpy.zeros(nsrc)
    #    for isrc in range(nsrc):
    #        samples[isrc]=u[y[isrc]-offset,x[isrc]-offset]
    #    sf='samples.txt'
    #    sf=os.path.join(dataset,sf)
    #    numpy.savetxt(sf,1.0e6*samples) # Big ~ 1 GB!
    #    backgroundS=numpy.median(samples)
    #    #backgroundS=0.0
        
    print 'Assuming a background mean flux = %f uJy' % backgroundS

    if injectRandomizeExtractedPosns:
        numpy.random.seed(seed=injectRandomizeExtractedPosnsSeed)
        x=numpy.random.randint(gaussianDim,imhead['NAXIS1']-gaussianDim,nsrc)
        y=numpy.random.randint(gaussianDim,imhead['NAXIS2']-gaussianDim,nsrc)


    for iS in range(nsrc):
        # Pixel photometry is cheap - but inaccurate
        if injectDoPosnNoise:
            g,mask=gaussian2(gaussianDim2,sig,rpix=rpix,\
                             xshift=xnoise[iS],yshift=ynoise[iS])
        Spix[iS]=Jy2muJy*(u[y[iS]+ynoise[iS],x[iS]+xnoise[iS]])-backgroundS
        pSpix[iS]=100.0*(Spix[iS]-S[iS])/S[iS]
        # Aperture photometry is slow -
        if injectAperturePhot:
            #SEx[iS]=Jy2muJy*psf_photometry(imdata[0,0,:,:],[(x[iS]+0.5,y[iS]+0.5)],PSF)[0]/fac/fwhm/fwhm
            apertureStamp=u[y[iS]-gaussianDim:y[iS]+gaussianDim+1,x[iS]-gaussianDim:x[iS]+gaussianDim+1]
            # g.sum() == fac*fwhm*fwhm to an excellent approximation
            ##SEx[iS]=2.0*Jy2muJy*numpy.sum(apertureStamp*g)/numpy.sum(g) # equivalent
            #SEx[iS]=2.0*Jy2muJy*numpy.average(apertureStamp,weights=g)#/fac/fwhm/fwhm
            ##SEx[iS]=Jy2muJy*numpy.mean(apertureStamp)*fac*fwhm*fwhm
            #print g.sum(),g[mask].sum()
            #SEx[iS]=Jy2muJy*numpy.mean(apertureStamp[mask])*fac*fwhm*fwhm
            ##SEx[iS]=Jy2muJy*numpy.sum(apertureStamp[mask]*g[mask])/g[mask].sum()
            ##SEx[iS]=Jy2muJy*numpy.average(apertureStamp[mask])*apercorr
            #Ssp[iS]=Jy2muJy*numpy.average(apertureStamp[mask],weights=g[mask])*g[mask].size*apercorr/g.sum() # g.sum() to normalize the gaussian
            #print g[mask].sum(), g.sum(),fac*fwhm*fwhm
            #print Jy2muJy*numpy.sum(apertureStamp*g)/numpy.sum(g)
            ##SEx[iS]=Jy2muJy*numpy.sum(apertureStamp[mask]*g[mask])/numpy.sum(g[mask])/erf2
            #Ssp[iS]=Jy2muJy*numpy.sum(apertureStamp[mask])/apertureStamp[mask].size
            #Ssp[iS]=2.65*Jy2muJy*apertureStamp[mask].sum()/apertureStamp[mask].size/erf2
            #Ssp[iS]=Jy2muJy*numpy.sum(apertureStamp[mask]*g[mask])/apertureStamp[mask].size/erf2
            #areaFac=pi*sig*sig
            ##areaFac=pi*rpix*rpix

            ##Ssp[iS]=Jy2muJy*numpy.sum(apertureStamp[mask])/erf2/areaFac
            ##SEx[iS]=Jy2muJy*numpy.average(apertureStamp[mask])/erf2/areaFac
            # gaussian-weighted sum inside aperture
            #u=upsample_image(apertureStamp,upsample_factor=4,use_numpy_fft=False)
            ##areaFac=pi*rpix*rpix
            ##SEx[iS]=Jy2muJy*numpy.average(apertureStamp[mask]*g[mask])*g[mask].size/erf2/(pi*sig*sig)/fac/0.9370906
            # This is just mean_<rpix{ai*g_i}*npix(<rpix) -> sum_<rpix{ai*g_i}
            # Then divide by erf2 -> total volume of gaussian for all pixels
            # Then normalize the gaussian
            # This = 0.986156348327 =~ 1 for rpix=0.5*u and sig=4*u/2.35482
            # The following areaFac gives SEx/Strue = 0.999598663473 ! :)
            # (for upsamplingFactor=8)
            areaFac=g[mask].size*(pi*rpix*rpix)/g[mask].sum()
            SEx[iS]=Jy2muJy*numpy.average(apertureStamp[mask]*g[mask])*(areaFac)/erf2/(2.0*pi*sig*sig)/fluxCorrFactor
            Smax[iS]=Jy2muJy*numpy.max(apertureStamp[mask])

            pSEx[iS]=100.0*(SEx[iS]-S[iS])/S[iS]
            #Ssp[iS]=Jy2muJy*numpy.average(apertureStamp[mask]*g[mask])/erf2
            print iS+1,rpix,sig,erf2,g[mask].sum(),g[mask].size,pi*rpix*rpix,g.size,Smax[iS],S[iS],Spix[iS],Spix[iS]/S[iS]
            line='%i %f %f %f %f %f %f %f %f %f %f %f %f'%(iS+1,x[iS],y[iS],S[iS],SEx[iS],SEx[iS]-S[iS],pSEx[iS],Spix[iS],Spix[iS]-S[iS],pSpix[iS],Ssp[iS],Smax[iS],numpy.mean(Spix[:iS]/S[:iS]))
#            print 'integral of gaussian mask = ',g[mask].sum()
#            print 'integral of gaussian stamp = ',g.sum()
#            print 'peak flux = ',S[iS]
#            print 'sum over stamp pixels = ',Jy2muJy*apertureStamp.sum()
#            print 'sum over masked pix = ',Jy2muJy*apertureStamp[mask].sum()
#            #print 'peak flux * integral = ',g[mask].sum()*S[iS]
#            print 'fac*fwhm*fwhm',fac*fwhm*fwhm
#            print 'pi*rpix*rpix',pi*rpix*rpix
#            print 'pi*sig*sig',pi*sig*sig
#            print 'Sum masked pix / pi*sig*sig',Jy2muJy*apertureStamp[mask].sum()/pi/sig/sig
#            print line
        ph.write('%s\n'%line)
#        sys.exit(0)
        if iS > 0 and iS % 100 == 0:
            print '\r==>',iS+1,x[iS],y[iS],S[iS],SEx[iS],SEx[iS]-S[iS],pSEx[iS],Spix[iS],pSpix[iS],#,numpy.mean(pS[:iS]),numpy.std(pS[:iS])
            ph.flush()
        #raw_input('')
    ph.close()
    print '\nPhotometry written to %s' % outphotf

    print 'SEx   %f +/- 1sigma=%f' % (numpy.mean(SEx[:nsrc]),numpy.std(SEx[:nsrc]))
    print 'Ssp   %f +/- 1sigma=%f' % (numpy.mean(Ssp[:nsrc]),numpy.std(Ssp[:nsrc]))
    print 'Spix  %f +/- 1sigma=%f' % (numpy.mean(Spix[:nsrc]),numpy.std(Spix[:nsrc]))
    print 'pSEx  %f +/- 1sigma=%f' % (numpy.mean(pSEx[:nsrc]),numpy.std(pSEx[:nsrc]))
    print 'pSpix %f +/- 1sigma=%f' % (numpy.mean(pSpix[:nsrc]),numpy.std(pSpix[:nsrc]))
    print 'bSEx %f +/- 1sigma=%f' % (numpy.mean(SEx[:nsrc]/S[:nsrc]),numpy.std(SEx[:nsrc]/S[:nsrc]))
    print 'bSpix %f +/- 1sigma=%f' % (numpy.mean(Spix[:nsrc]/S[:nsrc]),numpy.std(Spix[:nsrc]/S[:nsrc]))

    tail=SEx[numpy.where(SEx[:nsrc]<0.0)]
    symm=numpy.concatenate([tail,-tail])
    print 'SEx   %f +/- 1sigma=%f' % (numpy.mean(symm),numpy.std(symm))
    tail=Spix[numpy.where(Spix[:nsrc]<0.0)]
    symm=numpy.concatenate([tail,-tail])
    print 'Spix  %f +/- 1sigma=%f' % (numpy.mean(symm),numpy.std(symm))
    tail=Ssp[numpy.where(Ssp[:nsrc]<0.0)]
    symm=numpy.concatenate([tail,-tail])
    print 'Ssp   %f +/- 1sigma=%f' % (numpy.mean(symm),numpy.std(symm))

    # attempt to broadcast table to TOPCAT
    #from pylibs import mysamp
    #myclient=mysamp.demo()
    #x = numpy.arange(0, 100)
    #y=x**2
    #Stop and start TOPCAT connection to SAMP hub                                   
    #Broadcast table to TOPCAT                                                      
    #myclient['t0'] = {'x':x, 'y':y }
    #myclient.hub.stop()

    print 'Min Spix = %f uJy' % Spix.min()
    
    # dump map to file
    if not os.path.exists(injectedMapDir): os.mkdir(injectedMapDir)
    outmapf=os.path.join(injectedMapDir,injectedMap)
    pyfits.writeto(outmapf,u,header=imhead,clobber=True)
    print 'Injected map written to %s' % outmapf

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)


    
