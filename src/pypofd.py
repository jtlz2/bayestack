from scipy import *
import math
import numpy as np
from scipy import stats
from scipy import fftpack
import pylab as pyl

def pofd1(pxsz,logs,logdnds,mins,maxs,minx,maxx,bmaj,bmin,phi,sigma=0.,meanx=0.,summean=False):
    """ Calculated predicted P(D) given a source count and image parameters differs from pofd2 by calculating 2d elliptical gaussian for beam given input parameters. Assumes gaussian beam with no negative values or sidelobes.
    pxsz = pixel size in square arcseconds
    logs = log of flux in Jy
    logdnds = log of source count 
    mins = minimum s to use in Jy
    maxs =  max s to use in Jy
    minx = minimum histogram value in Jy/beam
    maxx = maximum histogram value in Jy/beam
    noise = instrumental noise value (same units as x) default is 0
    bmaj = beam major axis in arcseconds assuming gaussian beam
    bmin = beam minor axis in arcseconds assuming gaussian beam
    phi = beam position angle
    
    return : P(D), xx
    """

    pofd=0.
    pxsz_ns=sqrt(pxsz)
    assr=4.25452e10
    pii=np.pi
    #PIXEL SIZE IN STERADIANS
    pxsz_sr=pxsz/assr

    #LOG(S), DNDS, LOG(DNDS)
    ss=10.**logs
    #SCALE COUNTS TO PER PIXEL INSTEAD OF PER SR
    logdnds_px = logdnds+log10(pxsz_sr)

    #SETUP D ARRAY
    dnn = 2**18 #SIZE OF D ARRAY SOMETHING LARGE FOR GOOD RESOLUTION
    dx  =(maxx-minx)/(dnn-1.)
    dd = np.arange(minx,maxx+dx/2.,dx)
    dd2 = (dd[1:]+dd[:-1])/2.
    dnn1=len(dd)


    #SETUP BEAM ARRAY (ASSUME GAUSSIAN)

    bmaj_px = bmaj/pxsz_ns#BEAM FWHM IN PIXELS
    bmin_px = bmin/pxsz_ns#BEAM FWHM IN PIXELS
    bmajsd_px = bmaj_px/2.3548
    bminsd_px = bmin_px/2.3548
    phi_r = phi*pii/180.
    bx = np.arange(bmaj_px*3.) #CREATE X ARRAY TO CALC GAUSSIAN, PROBABLY 2 OR 3 TIMES THE SIZE
    bx = bx-(len(bx)-1)/2
    aa = (cos(phi_r)**2/(2.*bmajsd_px**2))+(sin(phi_r)**2/(2.*bminsd_px**2))
    bb = ((-1.)*sin(2*phi_r)**2/(4.*bmajsd_px**2))+(sin(2.*phi_r)**2/(4.*bminsd_px**2))
    cc = (sin(phi_r)**2/(2.*bmajsd_px**2))+(cos(phi_r)**2/(2.*bminsd_px**2))
    beam = exp((-1.)*(aa*bx.reshape(1,len(bx))**2+2.*bb*bx.reshape(1,len(bx))*bx.reshape(len(bx),1)+cc*bx.reshape(len(bx),1)**2))
    beam1 = beam.reshape(beam.size)
    beam1=beam1[beam1>=1.e-06] #SET SOME LOWER LIMIT
    beam11 = np.unique(beam1)
    nb=len(beam11)
    bcount=np.zeros(nb)
    for i in range(nb):
        bcount[i]=beam1[beam1==beam11[i]].size

    #DO RX INTEGRAL
    rx = np.zeros(dnn)
    xd = (np.arange(dnn)+1)*dx
    inti = ((log10(xd.max())-log10(dx))/500.)
    xdintl = np.arange(log10(dx),log10(xd.max()+inti/2.),inti)
    xdint=10.**xdintl
    ni=xdint.size
    rxint=np.zeros(ni)
    print mins,maxs
    for i in range(ni):
        flux = xdint[i]/beam11
        flux1 = flux[(flux>=mins)&(flux<maxs)]
        beam2 = beam11[(flux>=mins)&(flux<maxs)]
        bcnt2 = bcount[(flux>=mins)&(flux<maxs)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxint[i] = sum((10.**newy/beam2)*bcnt2)

    rx=np.interp(xd,xdint,rxint)

    #COMPUTE FFTS
    frx0=rx*dx
    frx1=fftpack.fft(frx0)
    frx2=frx1-frx1[0]

    omega=fftpack.fftfreq(dnn,dx)*np.pi*2

    ### IF Want to estimate mean ####
    if summean==True:
        summ=sum(xd*rx)
        mapmean=summ*dx
        mapmeani=complex(0,mapmean)
        frx2=frx2-(omega*mapmeani)

    ##IF YOU WANT TO ADD NOISE
    frx3=frx2-(sigma**2*omega**2)/2.

    frx4=exp(frx3)

    pofd=np.real(fftpack.ifft(frx4))

    cs=int(((meanx-dd[0])/dx))
    aa=np.argmax(pofd)
    sh=cs-aa
    if dd[0] !=0:
        pofd=np.roll(pofd,sh)
        
    
    return pofd,dd

def pofd2(pxsz,logs,logdnds,mins,maxs,minx,maxx,beam,sigma=0.,meanx=0.,summean=False):
    """ Calculated predicted P(D) given a source count and image parameters, differs from pofd1 by having user input array for beam, Assumes gaussian beam with no negative values or sidelobes.
    pxsz = pixel size in arcseconds
    logs = log of flux in Jy
    logdnds = log of source count 
    mins = minimum s to use in Jy
    maxs =  max s to use in Jy
    minx = minimum histogram value in Jy/beam
    maxx = maximum histogram value in Jy/beam
    noise = instrumental noise value (same units as x) default is 0
    beam = array containing values of beam peak normalized
    
    return : P(D), xx
    """
    pofd=0.


    assr=4.25452e10
    pii=np.pi
    #PIXEL SIZE IN STERADIANS
    pxsz_sr=pxsz/assr

    #LOG(S), DNDS, LOG(DNDS)
    ss=10.**logs
    #SCALE COUNTS TO PER PIXEL INSTEAD OF PER SR
    logdnds_px = logdnds+log10(pxsz_sr)

    #SETUP D ARRAY
    dnn = 2**18 #SIZE OF D ARRAY SOMETHING LARGE FOR GOOD RESOLUTION
    dx  =(maxx-minx)/(dnn-1.)
    dd = np.arange(minx,maxx+dx/2.,dx)
    dd2 = (dd[1:]+dd[:-1])/2.
    dnn1=len(dd)


    #SETUP BEAM ARRAY (ASSUME GAUSSIAN)
    jl=beam.shape
    jll=len(jl)
    if jll==2:
        beam1 = beam.reshape(beam.size)
    if jll==1:
        beam1 = beam.copy()
    beam1=beam1[beam1>=1.e-06] #SET SOME LOWER LIMIT
    beam11 = np.unique(beam1)
    nb=len(beam11)
    bcount=np.zeros(nb)
    for i in range(nb):
        bcount[i]=beam1[beam1==beam11[i]].size
    

    #DO RX INTEGRAL
    rx = np.zeros(dnn)
    xd = (np.arange(dnn)+1)*dx
    inti = ((log10(xd.max())-log10(dx))/500.)
    xdintl = np.arange(log10(dx),log10(xd.max()+inti/2.),inti)
    xdint=10.**xdintl
    ni=xdint.size
    rxint=np.zeros(ni)
    for i in range(ni):
        flux = xdint[i]/beam11
        flux1 = flux[(flux>=mins)&(flux<maxs)]
        beam2 = beam11[(flux>=mins)&(flux<maxs)]
        bcnt2 = bcount[(flux>=mins)&(flux<maxs)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxint[i] = sum((10.**newy/beam2)*bcnt2)

    rx=np.interp(xd,xdint,rxint)

    #COMPUTE FFTS
    frx0=rx*dx
    frx1=fftpack.fft(frx0)
    frx2=frx1-frx1[0]

    omega=fftpack.fftfreq(dnn,dx)*np.pi*2.

    ### IF Want to estimate mean ####
    if summean==True:
        summ=sum(xd*rx)
        mapmean=summ*dx
        mapmeani=complex(0,mapmean)
        frx2=frx2-(omega*mapmeani)

    ##IF YOU WANT TO ADD NOISE
    frx3=frx2-(sigma**2*omega**2)/2.

    frx4=exp(frx3)

    pofd=np.real(fftpack.ifft(frx4))

    cs=int(((meanx-dd[0])/dx))
    aa=np.argmax(pofd)
    sh=cs-aa
    if dd[0] !=0:
        pofd=np.roll(pofd,sh)
        
    
    return pofd,dd


def pofd22(pxsz,logs,logdnds,mins,maxs,minx,maxx,beam,bcount,sigma=0.,meanx=0.,summean=False):
    """ Calculated predicted P(D) given a source count and image parameters, differs from pofd1 by having user input array for beam with unique values and array of counts, Assumes gaussian beam with no negative values or sidelobes.
    pxsz = pixel size in arcseconds
    logs = log of flux in Jy
    logdnds = log of source count 
    mins = minimum s to use in Jy
    maxs =  max s to use in Jy
    minx = minimum histogram value in Jy/beam
    maxx = maximum histogram value in Jy/beam
    noise = instrumental noise value (same units as x) default is 0
    beam = array containing unique values of beam peak normalized
    bcount = array telling number of duplicates for beam
    
    return : P(D), xx
    """
    pofd=0.


    assr=4.25452e10
    pii=np.pi
    #PIXEL SIZE IN STERADIANS
    pxsz_sr=pxsz/assr

    #LOG(S), DNDS, LOG(DNDS)
    ss=10.**logs
    #SCALE COUNTS TO PER PIXEL INSTEAD OF PER SR
    logdnds_px = logdnds+log10(pxsz_sr)

    #SETUP D ARRAY
    dnn = 2**18 #SIZE OF D ARRAY SOMETHING LARGE FOR GOOD RESOLUTION
    dx  =(maxx-minx)/(dnn-1.)
    dd = np.arange(minx,maxx+dx/2.,dx)
    dd2 = (dd[1:]+dd[:-1])/2.
    dnn1=len(dd)


    #SETUP BEAM ARRAY (ASSUME GAUSSIAN)

    beam11=beam[beam>=1.e-06] #SET SOME LOWER LIMIT
    nb=len(beam11)


    #DO RX INTEGRAL
    rx = np.zeros(dnn)
    xd = (np.arange(dnn)+1)*dx
    inti = ((log10(xd.max())-log10(dx))/500.)
    xdintl = np.arange(log10(dx),log10(xd.max()+inti/2.),inti)
    xdint=10.**xdintl
    ni=xdint.size
    rxint=np.zeros(ni)
    for i in range(ni):
        flux = xdint[i]/beam11
        flux1 = flux[(flux>=mins)&(flux<maxs)]
        beam2 = beam11[(flux>=mins)&(flux<maxs)]
        bcnt2 = bcount[(flux>=mins)&(flux<maxs)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxint[i] = sum((10.**newy/beam2)*bcnt2)

    rx=np.interp(xd,xdint,rxint)

    #COMPUTE FFTS
    frx0=rx*dx
    frx1=fftpack.fft(frx0)
    frx2=frx1-frx1[0]

    omega=fftpack.fftfreq(dnn,dx)*np.pi*2.

    ### IF Want to estimate mean ####
    if summean==True:
        summ=sum(xd*rx)
        mapmean=summ*dx
        mapmeani=complex(0,mapmean)
        frx2=frx2-(omega*mapmeani)

    ##IF YOU WANT TO ADD NOISE
    frx3=frx2-(sigma**2*omega**2)/2.

    frx4=exp(frx3)

    pofd=np.real(fftpack.ifft(frx4))

    cs=int(((meanx-dd[0])/dx))
    aa=np.argmax(pofd)
    sh=cs-aa
    if dd[0] !=0:
        pofd=np.roll(pofd,sh)
        
    
    return pofd,dd

def pofd3(pxsz,logs,logdnds,mins,maxs,minx,maxx,beamc,beamd,cleanf,sigma=0.,meanx=0.,summean=False):
    """ Calculated predicted P(D) given a source count and image parameters, differs from pofd1 by having user input array for beam, differs from pofd1 and pofd2 by using dirty beam, with negative and positive sidelobes.
    pxsz = pixel size in arcseconds
    logs = log of flux in Jy
    logdnds = log of source count 
    mins = minimum s to use in Jy
    maxs =  max s to use in Jy
    minx = minimum histogram value in Jy/beam
    maxx = maximum histogram value in Jy/beam
    noise = instrumental noise value (same units as x) default is 0
    beamc = array containing values of clean beam peak normalized
    beamd = array containing values of dirty beam peak normalized
    cleanf = clean limit, flux below which dirty beam is used
    
    return : P(D), xx
    """
    pofd=0.
    

    assr=4.25452e10
    pii=np.pi
    #PIXEL SIZE IN STERADIANS
    pxsz_sr=pxsz/assr

    #LOG(S), DNDS, LOG(DNDS)
    ss=10.**logs
    #SCALE COUNTS TO PER PIXEL INSTEAD OF PER SR
    logdnds_px = logdnds+log10(pxsz_sr)

    #SETUP D ARRAY
    dnn = 2**18 #SIZE OF D ARRAY SOMETHING LARGE FOR GOOD RESOLUTION
    dx  =(maxx-minx)/(dnn-1.)
    dd = np.arange(minx,maxx+dx/2.,dx)
    dd2 = (dd[1:]+dd[:-1])/2.
    dnn1=len(dd)

    #SETUP BEAM ARRAY (ASSUME GAUSSIAN)
    jl=beamd.shape
    jll=len(jl)
    if jll==2:
        beam1d = beamd.reshape(beamd.size)
    if jll==1:
        beam1d = beamd.copy()
    jl=beamc.shape
    jll=len(jl)
    if jll==2:
        beam1c = beamc.reshape(beamc.size)
    if jll==1:
        beam1c = beamc.copy()
    beam1d=beam1d[abs(beam1d)>=1.e-06] #SET SOME LOWER LIMIT
    beam1c=beam1c[abs(beam1c)>=1.e-06] #SET SOME LOWER LIMIT
    beam11d = np.unique(beam1d)
    nbd=len(beam11d)
    bcountd=np.zeros(nbd)
    for i in range(nbd):
        bcountd[i]=beam1d[beam1d==beam11d[i]].size

    beam11c = np.unique(beam1c)
    nbc=len(beam11c)
    bcountc=np.zeros(nbc)
    for i in range(nbc):
        bcountc[i]=beam1c[beam1c==beam11c[i]].size
    
    bdpos=beam11d[beam11d>0]
    bdcpos=bcountd[beam11d>0]
    bdneg=beam11d[beam11d<0]
    bdcneg=bcountd[beam11d<0]

    nbdp=bdpos.size
    nbdn=bdneg.size
    #DO RX INTEGRAL
    rxc = np.zeros(dnn)
    rxdp = np.zeros(dnn)
    rxdn = np.zeros(dnn)

    xd = (np.arange(dnn)+1)*dx
    inti = ((log10(xd.max())-log10(dx))/500.)
    xdintl = np.arange(log10(dx),log10(xd.max()+inti/2.),inti)
    xdint=10.**xdintl
    ni=xdint.size
    rxintc=np.zeros(ni)
    rxintdp=np.zeros(ni)
    rxintdn=np.zeros(ni)
    #print ss.min(),ss.max(),cleanf
    #print beam11c.shape,bdpos.shape,bdneg.shape
    #print bdpos.min(),bdpos.max(),bdneg.min(),bdneg.max()
    #print bdpos.size()
    for i in range(ni):
        flux = xdint[i]/beam11c
        flux1 = flux[(flux>=mins)&(flux<maxs)]
        beam2 = beam11c[(flux>=mins)&(flux<maxs)]
        bcnt2 = bcountc[(flux>=mins)&(flux<maxs)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxintc[i] = sum((10.**newy/beam2)*bcnt2)

    for i in range(ni):
        flux = xdint[i]/bdpos
        flux1 = flux[(flux>=mins)&(flux<cleanf)]
        beam2 = bdpos[(flux>=mins)&(flux<cleanf)]
        bcnt2 = bdcpos[(flux>=mins)&(flux<cleanf)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxintdp[i] = sum((10.**newy/beam2)*bcnt2)

    for i in range(ni):
        flux = xdint[i]/abs(bdneg)
        flux1 = flux[(flux>=mins)&(flux<cleanf)]
        beam2 = abs(bdneg[(flux>=mins)&(flux<cleanf)])
        bcnt2 = bdcneg[(flux>=mins)&(flux<cleanf)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxintdn[i] = sum((10.**newy/beam2)*bcnt2)

    rxc=np.interp(xd,xdint,rxintc)
    rxdp=np.interp(xd,xdint,rxintdp)
    rxdn=np.interp(xd,xdint,rxintdn)

    #COMPUTE FFTS

    omega=fftpack.fftfreq(dnn,dx)*np.pi*2.

    frx0c=rxc*dx
    frx1c=fftpack.fft(frx0c)
    frx2c=frx1c-frx1c[0]

    frx0dp=rxdp*dx
    frx1dp=fftpack.fft(frx0dp)
    frx2dp=frx1dp-frx1dp[0]

    frx0dn=rxdn*dx
    frx1dn=fftpack.fft(frx0dn)
    frx2dn=frx1dn-frx1dn[0]
    n1p=np.arange(omega.size,dtype=int)
    os=np.argsort(omega)
    frx2dni=np.imag(frx2dn)
    frx2dni[os]=frx2dni[os]*-1.

    
    frx2=frx2c+frx2dp+frx2dn


    ### IF Want to estimate mean ####
    if summean==True:
        summ=sum(xd*rx)
        mapmean=summ*dx
        mapmeani=complex(0,mapmean)
        frx2=frx2-(omega*mapmeani)
    
    ##IF YOU WANT TO ADD NOISE
    frx3=frx2-(sigma**2*omega**2)/2.

    frx4=exp(frx3)

    pofd=np.real(fftpack.ifft(frx4))

    cs=int(((meanx-dd[0])/dx))
    aa=np.argmax(pofd)
    sh=cs-aa
    if dd[0] !=0:
        pofd=np.roll(pofd,sh)


    
    return pofd,dd


def pofd4(pxsz,logs,logdnds,mins,maxs,minx,maxx,beamc,bcountc,beamd,cntd,cleanf,sigma=0.,meanx=0.,summean=False):
    """ Calculated predicted P(D) given a source count and image parameters, differs from pofd1 by having user input array for beam, differs from pofd1 and pofd2 by using dirty beam, with negative and positive sidelobes and pofd3 by type of beam arrays
    pxsz = pixel size in arcseconds
    logs = log of flux in Jy
    logdnds = log of source count 
    mins = minimum s to use in Jy
    maxs =  max s to use in Jy
    minx = minimum histogram value in Jy/beam
    maxx = maximum histogram value in Jy/beam
    noise = instrumental noise value (same units as x) default is 0
    beamc = 1d array containing unique values of clean beam peak normalized
    cntc = 1d number of times to use each unique clean beam value
    beamd = 1d array containing unique values of dirty beam peak normalized
    cntd = 1d number of times to use each unique dirty beam value
    cleanf = clean limit, flux below which dirty beam is used
    
    return : P(D), xx
    """
    pofd=0.
    

    assr=4.25452e10
    pii=np.pi
    #PIXEL SIZE IN STERADIANS
    pxsz_sr=pxsz/assr

    #LOG(S), DNDS, LOG(DNDS)
    ss=10.**logs
    #SCALE COUNTS TO PER PIXEL INSTEAD OF PER SR
    logdnds_px = logdnds+log10(pxsz_sr)

    #SETUP D ARRAY
    dnn = 2**18 #SIZE OF D ARRAY SOMETHING LARGE FOR GOOD RESOLUTION
    dx  =(maxx-minx)/(dnn-1.)
    dd = np.arange(minx,maxx+dx/2.,dx)
    dd2 = (dd[1:]+dd[:-1])/2.
    dnn1=len(dd)

    #SETUP BEAM ARRAY (ASSUME GAUSSIAN)
    beam11d=beamd[abs(beamd)>=1.e-06] #SET SOME LOWER LIMIT
    beam11c=beamc[abs(beamc)>=1.e-06] #SET SOME LOWER LIMIT
    nbd=len(beam11d)
    nbc=len(beam11c)
    bdpos=beam11d[beam11d>0]
    bdcpos=bcountd[beam11d>0]
    bdneg=beam11d[beam11d<0]
    bdcneg=bcountd[beam11d<0]

    nbdp=bdpos.size
    nbdn=bdneg.size
    #DO RX INTEGRAL
    rxc = np.zeros(dnn)
    rxdp = np.zeros(dnn)
    rxdn = np.zeros(dnn)

    xd = (np.arange(dnn)+1)*dx
    inti = ((log10(xd.max())-log10(dx))/500.)
    xdintl = np.arange(log10(dx),log10(xd.max()+inti/2.),inti)
    xdint=10.**xdintl
    ni=xdint.size
    rxintc=np.zeros(ni)
    rxintdp=np.zeros(ni)
    rxintdn=np.zeros(ni)
    #print ss.min(),ss.max(),cleanf
    #print beam11c.shape,bdpos.shape,bdneg.shape
    #print bdpos.min(),bdpos.max(),bdneg.min(),bdneg.max()
    #print bdpos.size()
    for i in range(ni):
        flux = xdint[i]/beam11c
        flux1 = flux[(flux>=ss.min())&(flux<ss.max())]
        beam2 = beam11c[(flux>=ss.min())&(flux<ss.max())]
        bcnt2 = bcountc[(flux>=ss.min())&(flux<ss.max())]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxintc[i] = sum((10.**newy/beam2)*bcnt2)

    for i in range(ni):
        flux = xdint[i]/bdpos
        flux1 = flux[(flux>=ss.min())&(flux<cleanf)]
        beam2 = bdpos[(flux>=ss.min())&(flux<cleanf)]
        bcnt2 = bdcpos[(flux>=ss.min())&(flux<cleanf)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxintdp[i] = sum((10.**newy/beam2)*bcnt2)

    for i in range(ni):
        flux = xdint[i]/abs(bdneg)
        flux1 = flux[(flux>=ss.min())&(flux<cleanf)]
        beam2 = abs(bdneg[(flux>=ss.min())&(flux<cleanf)])
        bcnt2 = bdcneg[(flux>=ss.min())&(flux<cleanf)]
        logflx = log10(flux1)
        newy = np.interp(logflx,logs,logdnds_px)
        rxintdn[i] = sum((10.**newy/beam2)*bcnt2)

    rxc=np.interp(xd,xdint,rxintc)
    rxdp=np.interp(xd,xdint,rxintdp)
    rxdn=np.interp(xd,xdint,rxintdn)

    #COMPUTE FFTS

    omega=fftpack.fftfreq(dnn,dx)*np.pi*2.

    frx0c=rxc*dx
    frx1c=fftpack.fft(frx0c)
    frx2c=frx1c-frx1c[0]

    frx0dp=rxdp*dx
    frx1dp=fftpack.fft(frx0dp)
    frx2dp=frx1dp-frx1dp[0]

    frx0dn=rxdn*dx
    frx1dn=fftpack.fft(frx0dn)
    frx2dn=frx1dn-frx1dn[0]
    n1p=np.arange(omega.size,dtype=int)
    os=np.argsort(omega)
    frx2dni=np.imag(frx2dn)
    frx2dni[os]=frx2dni[os]*-1.

    
    frx2=frx2c+frx2dp+frx2dn


    ### IF Want to estimate mean ####
    if summean==True:
        summ=sum(xd*rx)
        mapmean=summ*dx
        mapmeani=complex(0,mapmean)
        frx2=frx2-(omega*mapmeani)
    
    ##IF YOU WANT TO ADD NOISE
    frx3=frx2-(sigma**2*omega**2)/2.

    frx4=exp(frx3)

    pofd=np.real(fftpack.ifft(frx4))

    cs=int(((meanx-dd[0])/dx))
    aa=np.argmax(pofd)
    sh=cs-aa
    if dd[0] !=0:
        pofd=np.roll(pofd,sh)


    
    return pofd,dd


def confnoise(xx,pofd1):
    """ Compute Confusion noise of given pofd

    returns confusion noise in microJy/beam (assuming input is in Jy/beam)
    """
    
    pofd_csum=np.cumsum(pofd1)
    ll=xx[np.argmin(abs(pofd_csum-.16))]
    uu=xx[np.argmin(abs(pofd_csum-.84))]
    conf=(uu-ll)/2.
    conf=conf/1.e-06

    return conf


def pdfint(xx,pofd1,newx):
    """ interpolate pofd to new axis
    xx = current pofd axis
    pofd = current pofd -- normalized such that sum(pofd*dx)=1
    newx = axis to interpolate tp
    
    returns newpofd
    """

    newpofd=np.interp(newx,xx,pofd1)
    

    return newpofd
