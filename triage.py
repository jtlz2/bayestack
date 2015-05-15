#!/usr/bin/env python

"""
"""

Jy2muJy=1.0e6
from math import e as euler
from math import log,log10,pi,sqrt
leuler=log10(euler)

import os,sys,numpy
import scipy
import pyfits
from scipy import sparse
from utils import bilinear_interpolation,sum2DGaussian
if os.getenv('PBS_O_HOST')==None or os.getenv('PBS_O_HOST').split('.')[0]!='baltasar':
    from photutils.psf import GaussianPSF
    from photutils import psf_photometry
    from skimage.morphology import disk,square,dilation

#-----------------------------------------------------------

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

def err(vector,alpha=0.68,version=2):

    """
    Calculate confidence intervals for a distribution
    DO NOT USE THIS FUNCTION FOR CALCULATING ERROR ON THE MEAN
    CURRENTLY SET UP SOLELY FOR THE MEDIAN (i.e. INCLUDES MAD FACTOR)
    """

    #from scipy.stats import morestats
    #import numpy
    #print numpy.median(vector)
    #(centre,(lim_low,lim_high))=\
    #    morestats.bayes_mvs(vector,alpha=alpha)[0]

    if version==1:
        # This is the (68 per cent) confidence interval about the median
        # Does not fall with sqrt(n)
        from scipy import stats

        percentile_median=50.0
        percentile_low=percentile_median-(100.0*alpha/2.0)
        percentile_high=percentile_median+(100.0*alpha/2.0)

        median   = stats.scoreatpercentile(vector,percentile_median)
        err_low  = median - stats.scoreatpercentile(vector,percentile_low)
        err_high = stats.scoreatpercentile(vector,percentile_high) - median
        #print err_low, err_high, median, percentile_low, percentile_high, percentile_median, (err_high+err_low)/2.0

    if version==2:
        # Take the error on the median to be MAD, i.e. sigma/1.4826
        # See e.g. http://en.wikipedia.org/wiki/Median_absolute_deviation#Relation_to_standard_deviation
        # NB This version assumes equally-weighted data!!
        import numpy
        from math import sqrt

        factor=1.4826

        sigma=numpy.std(vector)/sqrt(len(vector)-1)
        MAD=sigma/factor
        err_low=MAD
        err_high=MAD
        #print sigma,MAD,err_low,err_high

    return err_low,err_high,(err_high+err_low)/2.0

#-----------------------------------------------------------

def sumBivariateGaussian(sigma,npix):
    """
    Sum a bivariate gaussian, of width sigma, over npix x npix
    """

    import numpy
    from math import exp,pi,sqrt

    #sigma=1.7 # pix
    ic=jc=rc=(npix/2.0)-0.5

    x=numpy.zeros((npix,npix))

    for i in range(npix):
        for j in range(npix):
            prefactor=1.0/(2.0*pi*sigma*sigma)
            argument=-((float(i-ic))**2+(float(j-jc))**2)/(2.0*sigma*sigma)
            #print argument,prefactor
            x[i,j]=prefactor*exp(argument)
            #print i,j,x[i,j]

    peak_flux=numpy.max(x)

    # Normalize to unity peak
    x /= peak_flux
    print x

    sumx=numpy.sum(x)

    scaling_factor= 1.0/sumx

    return peak_flux,sumx,scaling_factor

#----------------------------------------------------------- 

def readCatalogue(do_shuffle=False,version=5,filename=None,cutKFlags=True,\
                  chibestthresh=50.0,cutHaloFlags=True,zReadThresh=5.0,\
                  reliabilityReadThresh=-100.0,IDRange=[0,-1]):

    """
    This reads the K-selected catalogue into memory (slow!)
    Cuts can be handled here if necessary.
    Version 4 is simply a counter of cut fractions
    Version 5 was added for bayestack.py project
    """

    import string,copy
    from random import shuffle

    #filename = \
    #    '/Users/jtlz2/video/matched-cats/full/head.10000'
#        '/Users/jtlz2/video/matched-cats/120203_nikhi/pbzk.txt'
    if version==1:
        filename = \
            '/Users/jtlz2/video/matched-cats/full/full_catalog_not-all-columns_with_flags.txt'
    elif version==2:
        filename = \
            '/Users/jtlz2/video/matched-cats/full-120504/VIDEO_04MAY_full.txt'
            #'/Users/jtlz2/video/matched-cats/full-120402/VIDEO_30MAR_full.txt'
            #'/Users/jtlz2/video/matched-cats/full-120307/VIDEO_07MAR_STELLARMASS_SFR.txt'
            #'/Users/jtlz2/video/matched-cats/full-120302/VIDEO_29FEB_STELLARMASS_SFR.txt'
            #'/Users/jtlz2/video/matched-cats/full-120223/VIDEO_23FEB_STELLARMASS_SFR.txt'
            #'/Users/jtlz2/video/matched-cats/full-with-masses-120220/VIDEO_20FEB_STELLARMASS.txt'
    #        '/Users/jtlz2/video/matched-cats/full-with-masses-120220/VIDEO_20FEB_STELLARMASS.txt'
    elif version==3 or version==4:
        filename = \
            '/Users/jtlz2/video/cats/120509/proc/VIDEO-120509.dat'

    elif version==5:
        filename=filename

# Format is different.. Identity, x, y, RA, dec, Kmag, Kmag_error,
# Kflag, haloflag, chi_QSO, chi_star, zbest, zQSO, zbest_68low,
# zbest_68high


    print 'Reading %s (version %i)...' % (filename,version)

    file=open(filename)
    clist=file.readlines()
    file.close()

    clist = [line[:-1] for line in clist]
    # Delete header
    if string.count(clist[0],'#') > 0:
        del clist[0]

    gals=[]

    if version==1:

        print 'version %i no longer supported --- aborting' % version
        return

        for n,line in enumerate(clist):
            # Split on whitespace
            l=line.split()

            reliability=float(l[18]) # < 0 non-detections > 0 detections
            z=float(l[13]) # z < 0.0 error and z >= 1000.0 error

            if reliability < 1.0 and z > 0.0 and z < 1000.0:
                gal=galaxy(l[0],l[1],l[2],l[3],l[4],l[5],l[6],l[7],l[8],l[9],\
                           l[10],l[11],l[12],l[13],l[14],l[15],l[16],l[17],l[18],\
                           l[19],l[20],l[21],l[22],l[23],filename,version=version)
                gals.append(gal)
            if n % 5000 == 0: print '==> %6i <==' % n

    elif version==2:

        chi_best_thresh=50.0
        k_class_star_thresh=1.10
        z_thresh=1000.0
        bmass_thresh=1.0e6
        #gmass_thresh=-90.0

        for n,line in enumerate(clist):
            # Split on whitespace
            l=line.split()

            z=float(l[17])             # z < 0.0 error and z >= 1000.0 error
            chibest=float(l[14])       # SED poor fit
            k_flags=int(float(l[11]))   # 0 = OK >0 not
            k_class_star=float(l[12]) # Assume > 0.90 is a star or QSO
            halo=int(float(l[13]))      # 1 = halo 0 = not
            #reliability=float(l[22])   # < 0 non-detections > 0 detections -- NOW CUT LATER
            #gmass=float(l[41])         # Mass from Elodie (flag -90)
            bmass=float(l[35])          # Mass from Dave (flag 1.0e6/7)


            if z > 0.0 and z < z_thresh \
                and halo==0 and k_flags==0 \
                and chibest < chi_best_thresh \
                and bmass < bmass_thresh and k_class_star < k_class_star_thresh:
                #and gmass > gmass_thresh \
                gal=galaxy(l,filename=filename,version=version)
                gals.append(gal)
            if n % 5000 == 0: print '==> %6i <==' % n


    elif version==3:

        chi_best_thresh=50.0
        k_class_star_thresh=1.10
        z_thresh=5.0
        bmass_thresh=1.0e6
        #gmass_thresh=-90.0

        for n,line in enumerate(clist):
            # Split on whitespace
            l=line.split()

            z=float(l[22])             # z < 0.0 error and z >= 1000.0 error
            chibest=float(l[19])       # SED poor fit
            k_flags=int(float(l[16]))   # 0 = OK >0 not
            k_class_star=float(l[17]) # Assume > 0.90 is a star or QSO
            halo=int(float(l[18]))      # 1 = halo 0 = not
            #reliability=float(l[--])   # < 0 non-detections > 0 detections -- NOW CUT LATER
            #gmass=float(l[--])         # Mass from Elodie (flag -90)
            bmass=float(l[31])          # Mass from Dave (flag 1.0e6/7)

            if z > 0.0 and z < z_thresh \
                and halo==0 and k_flags==0 \
                and chibest < chi_best_thresh \
                and bmass < bmass_thresh and k_class_star < k_class_star_thresh:
                #and gmass > gmass_thresh \
                gal=galaxy(l,filename=filename,version=version)
                gals.append(gal)
            if n % 5000 == 0: print '==> %6i <==' % n

        #if True or l[19]>0.0:
        #    gals.append(gal)

    elif version==5:

        chi_best_thresh=chibestthresh
        k_class_star_thresh=1.10
        z_thresh=1.0e10
        bmass_thresh=1.0e100 # was 1.0e6
        #gmass_thresh=-90.0
        rel_thresh=reliabilityReadThresh
        
        for n,line in enumerate(clist):
            # Split on whitespace
            l=line.split()
            #print line
            ident=int(l[0])           # ID number
            z=float(l[10])             # z < 0.0 error and z >= 1000.0 error
            chibest=float(l[7])       # SED poor fit
            k_flags=int(float(l[25]))   # 0 = OK >0 not
            k_class_star=float(l[26]) # Assume > 0.90 is a star or QSO
            halo=int(float(l[27]))      # 1 = halo 0 = not
            reliability=float(l[37])   # < 0 non-detections > 0 detections
            #gmass=float(l[--])         # Mass from Elodie (flag -90)
            #JK=float(l[19])
            #gi=float(l[20])
            try:
                bmass=float(l[15])          # Mass from Dave (flag 1.0e6/7)
            except ValueError:
                # This should no longer be needed since missing values
                # have been replaced via radiomatch.stil
                print 'Dummied a missing mass value'
                bmass=bmass_thresh+1.0

            if not cutKFlags: k_flags=0
            if not cutHaloFlags: halo=0
            ztest = (z > 0.0 and z < z_thresh)
            ztest=True
            if ztest \
                and ident < IDRange[-1] \
                and halo==0 and k_flags==0 \
                and chibest < chi_best_thresh \
                and bmass < bmass_thresh \
                and k_class_star < k_class_star_thresh \
                and reliability >= rel_thresh:
                #and gmass > gmass_thresh \
                gal=galaxy(l,filename=filename,version=version)
                gals.append(gal)
            if n % 5000 == 0: print '\r==> %6i <==' % n,

    elif version==4:

        print "This is a cut-fraction counter"

        chi_best_thresh=50.0
        k_class_star_thresh=1.10
        z_thresh=zReadThresh # Was 5.0
        bmass_thresh=1.0e6
        #gmass_thresh=-90.0

        n=0; nkeep=0; nz=0; nchi=0; nb=0; nhalo=0; nkf=0;
        for n,line in enumerate(clist):
            # Split on whitespace
            l=line.split()

            z=float(l[22])             # z < 0.0 error and z >= 1000.0 error
            chibest=float(l[19])       # SED poor fit
            k_flags=int(float(l[16]))   # 0 = OK >0 not
            #k_class_star=float(l[17]) # Assume > 0.90 is a star or QSO
            halo=int(float(l[18]))      # 1 = halo 0 = not
            #reliability=float(l[36])   # < 0 non-detections > 0 detections -- NOW CUT LATER
            #gmass=float(l[--])         # Mass from Elodie (flag -90)
            bmass=float(l[31])          # Mass from Dave (flag 1.0e6/7)

            if z < 0.0 or z > z_thresh: nz+=1
            if chibest > chi_best_thresh: nchi+=1
            if halo != 0: nhalo+=1
            if k_flags!=0: nkf+=1
            if bmass > bmass_thresh: nb+=1

            if z > 0.0 and z < z_thresh \
                and halo==0 and k_flags==0 \
                and chibest < chi_best_thresh \
                and bmass < bmass_thresh:
                nkeep+=1

            if n % 5000 == 0: print '==> %6i <==' % n

        print n,nkeep,nz,nchi,nhalo,nkf,nb
        # z<3:
        # import numpy
        #y=numpy.array([429806.0,278179.0,60110.0,36397.0,16720.0,69857.0,459.0])
        #100.0*y/y[0]        

        return None

    #def __init__(self,id,ximage,yimage,alpha,delta,kmag,errkmag,kflag,\
    #             kclass_star,haloflag,chibest,chiqso,chistar,z_best,\
    #             z_qso,z_low,z_high,flag,flag_pbzk,flag_sbzk,\
    #             flag_star,flag_ero,flag_drg,flag_normal):

    print 'Read %i objects' % len(gals)

    if do_shuffle:
        shuffle(gals)
        print 'Entries shuffled'
    else:
        print 'Entries not shuffled'

    return gals

#-----------------------------------------------------------

def selectGals(gals,flag='all',z_cut_low=0.0,z_cut_high=5.0,\
               k_cut_faint=23.5,k_cut_bright=0.0,\
               M_stellar_low=-1.0e10,M_stellar_high=1.0e10,\
               kabs_cut_faint=0.0,kabs_cut_bright=-100.0,masses='gio',\
               include_detections=False):
    """
    Make various inline cuts to the galaxy catalogue
    Some cuts (haloflag, SEx K_flags, reliability, poor z fits) \
        are handled at the read level (to cut down on the number of sources)
    """

    print 'selecting %s' % flag

    # non-detections
    # Currently handled in readCatalogue..
    #gal.flag<0.0
    # detections
    # Currently handled in readCatalogue..
    #gal.flag>0.0


    # ======== ABANDONED ========
    # Select sub-sample species
    if flag=='sbzk':
        gg = [gal for gal in gals if gal.flag_sbzk == 1]

    elif flag=='pbzk':
        gg = [gal for gal in gals if gal.flag_pbzk == 1]

    elif flag=='nbzk':
        gg = [gal for gal in gals if (gal.flag_sbzk < 1 and gal.flag_pbzk < 1 \
                                       and gal.flag_star_dunne < 1)]

    elif flag=='sgzk':
        gg = [gal for gal in gals if gal.flag_nsgzk == 1]

    elif flag=='pgzk':
        gg = [gal for gal in gals if gal.flag_npgzk == 1]

    elif flag=='ngzk':
        gg = [gal for gal in gals if (gal.flag_nsgzk < 1 and gal.flag_npgzk < 1 \
                                      #and gal.flag_star_dunne < 1)]
                                       and gal.flag_locus < 1)]


    elif flag=='sbc':
        gg = [gal for gal in gals if (gal.mod_best > 21 and gal.mod_best <= 37\
                                      and gal.flag_locus < 1)]

    elif flag=='scd':
        gg = [gal for gal in gals if (gal.mod_best > 37 and gal.mod_best <= 48\
                                      and gal.flag_locus < 1)]

    elif flag=='irr':
        gg = [gal for gal in gals if (gal.mod_best > 48 and gal.mod_best <= 58\
                                      and gal.flag_locus < 1)]

    elif flag=='eros':
        gg = [gal for gal in gals if gal.flag_ero == 1]

    elif flag=='drgs':
        #print [gal.flag_drg for gal in gals if gal.flag_drg==1]
        gg = [gal for gal in gals if gal.flag_drg == 1]
    # ===========================


    # Ellipticals
    elif flag=='ell':
        gg = [gal for gal in gals if (gal.flag_ell > 0 and gal.flag_locus_baldry > 0)]
    # mid (Sbc, Sdc, Im)
    elif flag=='mid':
        gg = [gal for gal in gals if (gal.flag_mid > 0 and gal.flag_locus_baldry > 0)]
    # Starburst galaxies
    elif flag=='sbn':
        gg = [gal for gal in gals if (gal.flag_sbn > 0 and gal.flag_locus_baldry > 0)]
    # These two are synonyms:
    elif flag=='not_star' or flag=='all':
        gg = [gal for gal in gals if gal.flag_locus_baldry > 0]       # Baldry

    elif flag=='star':
        gg = [gal for gal in gals if gal.flag_locus_baldry < 1]       # Baldry

    else:
        print 'unrecognized cut %s' % flag
        flag=='unrecon'
        gg=gals

    #print 'gg',len(gg),M_stellar_low,M_stellar_high,z_cut_low,z_cut_high
    # Apply cuts on z_best, kmag_aper_2, mass_best and kabsmag
    print 'Keeping %4.2f < z < %4.2f + %4.2f < kmag_aper_2 < %4.2f + %4.2f < mass_best < %4.2f + %4.2f < kabsmag < %4.2f' \
        % (z_cut_low,z_cut_high,k_cut_bright,k_cut_faint,\
           M_stellar_low,M_stellar_high,kabs_cut_bright,kabs_cut_faint)

    # Eliminate objects with any of g, z or K missing magnitudes
    #gh = [gal for gal in gg if gal.flag_normal > -1]
    #print 'excising -1 flags'
    gh=gg

    # Eliminate McA detections, if required
    #if include_detections:
    #    print 'including McA sources'
    #    gi=gh
    #else:
    #    print 'excluding McA sources'
    #    gi=[gal for gal in gh if gal.flag < 0]
    gi=gh

    # Make the other cuts
    #masses='bon'
    if masses=='gio':
        gj = [gal for gal in gi if (gal.z_best<z_cut_high and gal.z_best>=z_cut_low \
                                    and gal.kmag<k_cut_faint and gal.kmag>=k_cut_bright \
                                    and gal.mass_gio<M_stellar_high \
                                    and gal.mass_gio>=M_stellar_low) \
                                    and gal.kabsmag<kabs_cut_faint \
                                    and gal.kabsmag>=kabs_cut_bright ]
    elif masses=='bon':
        gj = [gal for gal in gi if (gal.z_best<z_cut_high and gal.z_best>=z_cut_low \
                                    and gal.kmag<k_cut_faint and gal.kmag>=k_cut_bright \
                                    and gal.mass_bon<M_stellar_high \
                                    and gal.mass_bon>=M_stellar_low) \
                                    and gal.kabsmag<kabs_cut_faint \
                                    and gal.kabsmag>=kabs_cut_bright ]

    else:
        print 'invalid mass type! --- %s' % masses
        return

    print 'mass type %s' % masses

    # z < 0.2  K <= 23.0
    #gg = [gal for gal in gals if (gal.z_best<0.2 and gal.kmag<=23.0)]
    # z < 0.7  K <= 23.0
    #gg = [gal for gal in gals if (gal.z_best<0.7 and gal.kmag<=23.0)]
    # z > 0.7  K <= 23.0
    #gg = [gal for gal in gals if (gal.z_best>0.7 and gal.kmag<=23.0)]
    # z > 0.2  K <= 21.5
    #gg = [gal for gal in gals if (gal.z_best>0.2 and gal.kmag<=21.5)]
    # z > 0.2  21.5 < K <= 23.0
    #gg = [gal for gal in gals if (gal.z_best>0.2 and gal.kmag>21.5 and gal.kmag<=23.0)]

    del gg,gh,gi

    return gj

#----------------------------------------------------------- 

class galaxy(object):
    """
    This holds the properties of each K-selected source
    """

    def __init__(self,l,\
    #(id,ximage,yimage,alpha,delta,kmag,errkmag,kflag,
    #            kclass_star,haloflag,chibest,chiqso,chistar,z_best,
    #            z_qso,z_low,z_high,mod_best,flag,flag_pbzk,flag_sbzk,
    #            flag_star,flag_star1,flag_qso,flag_locus,flag_star_dunne,
    #            flag_ero,flag_drg,flag_normal,bonf_massbest,bonf_masserr,
    #            bonf_sfrbest,bonf_sfrerr,gio_sfr,gio_sfrerr,gio_mass,
    #            gio_masserr),
                 filename=None,version=5):

        # Added for bayestack.py work
        if version==5:
            self.id=int(l[0])
            self.alpha=float(l[1])                   # 3
            self.delta=float(l[2])                   # 4
            self.kmag=float(l[3])                   # 15
            #self.kmag=float(l[28])                   # [K_MAG_AUTO]
            self.kflag=int(float(l[25]))             # 16
            self.kclass_star=float(l[26])            # 17
            self.haloflag=int(float(l[27]))          # 18
            self.chibest=float(l[7])                # 19
            self.z_best=float(l[10])                 # 22
            self.mod_best=int(float(l[14]))           # 29
            self.kabsmag=float(l[18])                # 30

            self.mass_bon=float(l[15])               # 31
            self.err_mass_bon=-99.0           # 32
            self.sfr_bon=-99.0                # 35
            self.err_sfr_bon=-99.0            # 36
            self.flag_ell=int(float(l[21]))          # 37
            self.flag_mid=int(float(l[22]))          # 38
            self.flag_sbn=int(float(l[23]))          # 39

            self.flag_locus_baldry=int(float(l[24])) # 42

            self.filename=filename                   # 43
            self.lumdist=0.0                         # Placeholder
            self.factor=0.0                          # Placeholder
            self.age=0.0                             # Placeholder

            self.sfr_gio=1.0                         # Dummy
            self.err_sfr_gio=1.0                     # Dummy
            self.mass_gio=1.0                        # Dummy
            self.err_mass_gio=1.0                    # Dummy

            self.radio_pflux=float(l[30])            # Radio properties
            self.radio_tflux=float(l[32])
            self.radio_sep=float(l[34])
            self.reliability=float(l[37])
            self.alpha_radio=float(l[38])
            self.delta_radio=float(l[39])

            self.flag=self.reliability
            self.flag_pbzk=1
            self.flag_sbzk=1
            self.flag_npgzk=1
            self.flag_nsgzk=1
            self.flag_star=1
            self.flag_star1=1
            self.flag_qso=1
            self.flag_locus=1
            self.flag_star_dunne=1
            self.flag_ero=1
            self.flag_drg=1
            self.flag_normal=1
            #self.pflux=float(l[29])
            #self.pfluxerr=float(l[30])
            #self.tflux=float(l[31])
            #self.tfluxerr=float(l[32])


        elif version==3:
            ##                                       # l[?]
            self.id=int(l[0])                        # 0
            #self.ximage=float(l[1])                  # 1
            #self.yimage=float(l[2])                  # 2
            self.alpha=float(l[3])                   # 3
            self.delta=float(l[4])                   # 4
            #self.umag=float(l[5])                    # 5
            #self.gmag=float(l[6])                    # 6
            #self.rmag=float(l[7])                    # 7
            #self.imag=float(l[8])                    # 8
            #self.i2mag=float(l[9])                   # 9
            #self.zmag=float(l[10])                   # 10
            #self.z2mag=float(l[11])                  # 11
            #self.ymag=float(l[12])                   # 12
            #self.jmag=float(l[13])                   # 13
            #self.hmag=float(l[14])                   # 14
            self.kmag=float(l[15])                   # 15
            self.kflag=int(float(l[16]))             # 16
            self.kclass_star=float(l[17])            # 17
            self.haloflag=int(float(l[18]))          # 18
            self.chibest=float(l[19])                # 19
            #self.chiqso=float(l[20])                 # 20
            #self.chistar=float(l[21])                # 21
            self.z_best=float(l[22])                 # 22
            #self.z_qso=float(l[23])                  # 23
            #self.z_low=float(l[24])                  # 24
            #self.z_high=float(l[25])                 # 25

            #self.z_ml=float(l[26])                   # 26
            #self.z_ml_low=float(l[27])               # 27
            #self.z_ml_high=float(l[28])              # 28
            self.mod_best=int(float(l[29]))           # 29
            self.kabsmag=float(l[30])                # 30

            self.mass_bon=float(l[31])               # 31
            self.err_mass_bon=float(l[32])           # 32
            #self.z_best_2=float(l[33])               # 33
            #self.z_best_err_2=float(l[34])           # 34
            self.sfr_bon=float(l[35])                # 35
            self.err_sfr_bon=float(l[36])            # 36
            self.flag_ell=int(float(l[37]))          # 37
            self.flag_mid=int(float(l[38]))          # 38
            self.flag_sbn=int(float(l[39]))          # 39
            #self.gi=float(l[40])                     # 40
            #self.jk=float(l[41])                     # 41
            self.flag_locus_baldry=int(float(l[42])) # 42

            self.filename=filename                   # 43
            self.lumdist=0.0                         # Placeholder
            self.factor=0.0                          # Placeholder
            self.age=0.0                             # Placeholder

            self.sfr_gio=1.0                         # Dummy
            self.err_sfr_gio=1.0                     # Dummy
            self.mass_gio=1.0                        # Dummy
            self.err_mass_gio=1.0                    # Dummy
            self.flag=1
            self.flag_pbzk=1
            self.flag_sbzk=1
            self.flag_npgzk=1
            self.flag_nsgzk=1
            self.flag_star=1
            self.flag_star1=1
            self.flag_qso=1
            self.flag_locus=1
            self.flag_star_dunne=1
            self.flag_ero=1
            self.flag_drg=1
            self.flag_normal=1


        if version==2:
            ##                                      # l[?]
            self.id=int(l[0])                       # 0
            #self.ximage=float(l[1])                 # 1
            #self.yimage=float(l[2])                 # 2
            self.alpha=float(l[3])                  # 3
            self.delta=float(l[4])                  # 4
            self.kmag=float(l[5])                   # 5
            #self.errkmag=float(l[6])                # 6
            #self.gmag=float(l[7])                   # 7
            #self.errgmag=float(l[8])                # 8
            #self.zmag=float(l[9])                   # 9
            #self.errzmag=float(l[10])               # 10

            self.kflag=float(l[11])                  # 11
            self.kclass_star=float(l[12])            # 12
            self.haloflag=float(l[13])               # 13

            self.chibest=float(l[14])               # 14
            self.chiqso=float(l[15])                # 15
            self.chistar=float(l[16])               # 16
            self.z_best=float(l[17])                # 17
            self.z_qso=float(l[18])                 # 18
            self.z_low=float(l[19])                 # 19
            self.z_high=float(l[20])                # 20
            self.mod_best=float(l[21])              # 21

            self.flag=float(l[22])                  # 22
            self.flag_pbzk=int(float(l[23]))        # 23
            self.flag_sbzk=int(float(l[24]))        # 24
            self.flag_npgzk=int(float(l[25]))       # 25
            self.flag_nsgzk=int(float(l[26]))       # 26
            self.flag_star=int(float(l[27]))        # 27
            self.flag_star1=int(float(l[28]))       # 28    # What is this!?
            self.flag_qso=int(float(l[29]))         # 29
            self.flag_locus=int(float(l[30]))       # 30
            self.flag_star_dunne=int(float(l[31]))  # 31
            self.flag_ero=int(float(l[32]))         # 32
            self.flag_drg=int(float(l[33]))         # 33
            self.flag_normal=int(float(l[34]))      # 34

            self.mass_bon=float(l[35])              # 35
            self.err_mass_bon=float(l[36])          # 36
            self.sfr_bon=float(l[37])               # 37
            self.err_sfr_bon=float(l[38])           # 38

            self.sfr_gio=float(l[39])               # 39
            self.err_sfr_gio=float(l[40])           # 40
            self.mass_gio=float(l[41])              # 41
            self.err_mass_gio=float(l[42])          # 42
            #self.sfr10_gio=float(l[43])             # 43
            #self.err_sfr10_gio=float(l[44])         # 44
            #self.srf100_gio=float(l[45])            # 45
            #self.err_sfr100_gio=float(l[46])        # 46
            #self.sfr_X2_gio=float(l[47])            # 47
            #self.sfr10_X2_gio=float(l[48])          # 48
            #self.sfr100_X2_gio=float(l[49])         # 49
            self.massX_gio=float(l[50])             # 50
            self.X_gio=float(l[51])                 # 51 Chi^2
            self.kabsmag=float(l[52])               # 52
            #self.jmag=float(l[53])                  # 53
            #self.errjmag=float(l[54])               # 54

            self.filename=filename                   # 35
            self.lumdist=0.0                         # Placeholder
            self.factor=0.0                          # Placeholder
            self.age=0.0                             # Placeholder

        if version==1:
            ##                                       # l[?]
            self.id=int(id)                          # 0
            #self.ximage=float(ximage)               # 1
            #self.yimage=float(yimage)               # 2
            self.alpha=float(alpha)                  # 3
            self.delta=float(delta)                  # 4
            self.kmag=float(kmag)                    # 5
            #self.errkmag=float(errkmag)             # 6
            self.kflag=float(kflag)                  # 7
            self.kclass_star=float(kclass_star)      # 8
            self.haloflag=float(haloflag)            # 9
            self.chibest=float(chibest)              # 10
            self.chiqso=float(chiqso)                # 11
            self.chistar=float(chistar)              # 12
            self.z_best=float(z_best)                # 13
            self.z_qso=float(z_qso)                  # 14
            self.z_low=float(z_low)                  # 15
            self.z_high=float(z_high)                # 16
            self.flag=float(flag)                    # 17
            self.flag_pbzk=int(float(flag_pbzk))     # 18
            self.flag_sbzk=int(float(flag_sbzk))     # 19
            self.flag_star=int(float(flag_star))     # 20
            self.flag_ero=int(float(flag_ero))       # 21
            self.flag_drg=int(float(flag_drg))       # 22
            self.flag_normal=int(float(flag_normal)) # 23
            self.filename=filename                   # 24
            self.lumdist=0.0                         # Placeholder
            self.factor=0.0                          # Placeholder

        return

#-----------------------------------------------------------

def running_avg(nums):
    """
    See http://stackoverflow.com/questions/1790550/running-average-in-python
    """
    for count in xrange(1, len(nums)+1):
        if count % 100 == 0: print '==> %6i <==' % count
        yield sum(nums[:count])/count

    return

#-----------------------------------------------------------

def online_variance(datum,n=0.0,mean=0.0,M2=0.0):
    """
    See http://stackoverflow.com/questions/5543651/computing-standard-deviation-in-a-stream
    """
    n += 1
    delta = datum - mean
    mean = mean + delta/n
    M2 = M2 + delta*(datum - mean)  
    variance = M2/n    # use this for population variance
    # variance = M2/(n - 1)  # use this for sample variance
    return variance,n,mean,M2

#-----------------------------------------------------------

class cube():

    """
    """

    def __init__(self,gals,radio_map,radio_header,simulate):
        """
        Generates the individual S and L cutouts
        Data simulation is controlled from here if needed
        """

        import cosmocalc
        from math import pi
        import copy
        import numpy
        from numpy.random import RandomState

        self.gals=gals
        del gals

        self.simulate=simulate

        if self.simulate:
            seed=1234
            prng=RandomState(seed)
        else:
            prng=None

        self.nx=self.ny=41 # 41

        self.stamps=numpy.zeros((numpy.size(self.gals),self.nx,self.ny))
        self.stamps_L=numpy.zeros((numpy.size(self.gals),self.nx,self.ny))
        print numpy.shape(self.stamps)

        for n,gal in enumerate(self.gals):
            st=stamp(gal,radio_map,radio_header,self.nx,self.ny,self.simulate,prng)
            self.stamps[n,:,:]=st.cutout

            # Now calculate L map = f(z)
            lumdist = cosmocalc.cosmocalc(gal.z_best)['DL_Mpc']
            factor = 9.52e18 * lumdist * lumdist * 4 * pi * (1+gal.z_best)**(-0.2)
            L=copy.deepcopy(st)
            L.cutout = factor * st.cutout

            self.stamps_L[n,:,:]=L.cutout
            del st, L, factor, lumdist
            if n % 200 == 0: print '==> %6i <==' % n

        return

#-----------------------------------------------------------

class sprint(object):
    """
    """

    def __init__(self,gals,seed=None):
        """
        Don't do very much - so can reload quickly
        Set up the random seed [actually even if not required]
        This is just a container class for gals
        """
        from numpy.random import RandomState

        self.gals=gals
        del gals

        if seed is None:
            self.seed=1234
        else:
            self.seed=seed
        print '[seed is %i] ---' % self.seed,
        self.prng=RandomState(self.seed)
        return


    def makeFrames(self,radio_map,radio_header,clip=True):
        """
        Generate the frame images and the mean stack (Dunne Fig. 2)
        """

        import numpy
        import matplotlib.pyplot as plt

        #self.gals=gals
        #del gals

        self.nx=self.ny=41

        clipn=5
        m=2000

        # Need to add the noise-weighted mean

        f=open('raw.txt','w')
        f.write('# n mean_S median_S max_S std_S mean_L median_L max_L std_L\n')

        self.simulate=False; prng=None
        #x=0.0; y=0.0; z=0.0
        nn=0
        run_std_S=[]
     #   means=[]; medians=[]
        frame_S=numpy.zeros((self.nx,self.ny))
        frame_L=numpy.zeros((self.nx,self.ny))
        for n,gal in enumerate(self.gals):
            # Extract cutout for this galaxy
            self.stampp=stamp(gal,radio_map,radio_header,self.nx,self.ny,self.simulate,prng)

            # Add the cutout to the summed frame, but only if not clipped
            if clip and (self.stampp.sovern < clipn and self.stampp.sovern > -clipn):
                frame_S += self.stampp.cutout
                frame_L += self.stampp.cutout_L
                nn+=1

            # Write results every m frames
            if n > 0 and n % m == 0:
                print n,
                mean_frame_S=frame_S/float(nn)
                self.writeFrame(mean_frame_S,nn,flux='S')
                print n,
                mean_frame_L=frame_L/float(nn)
                self.writeFrame(mean_frame_L,nn,flux='L')
                # Calculate the mean, median and std of the cutout
                run_std_S.append(numpy.std(mean_frame_S))
                run='%i %f %f %f %f %e %e %e %e\n' \
                    % (n,numpy.average(mean_frame_S),numpy.median(mean_frame_S),
                       numpy.max(mean_frame_S),numpy.std(mean_frame_S),\
                       numpy.average(mean_frame_L),numpy.median(mean_frame_L),
                       numpy.max(mean_frame_L),numpy.std(mean_frame_L))
                f.write(run)
                f.flush()


        f.close()

        # Now write the final mean frame
        frame_S /= float(nn)
        frame_L /= float(nn)
        self.writeFrame(frame_S,nn,f='mean_S.png',flux='S')
        self.writeFrame(frame_L,nn,f='mean_L.png',flux='L')

        self.writeSTD(run_std_S,m)

            # Care here - what is n if not clipped??
       #     else:
       #         nn=n
#            variance,n,mean,M2 = online_variance(self.stampp.mean,n,mean,M2)
#            online_std=Jy2muJy*numpy.sqrt(variance)
#            online_mean=Jy2muJy*mean
#            median=self.stampp.median
#            run='%i %f %f %f\n' % (n,online_mean,Jy2muJy*median,online_std)
#            print run
            #running_means.append(online_mean)
            #running_stds.append(online_std)

#            print Jy2muJy*self.stampp.mean,Jy2muJy*self.stampp.median
#            means.append(Jy2muJy*self.stampp.mean)
#            medians.append(Jy2muJy*self.stampp.median)
         #   x += Jy2muJy*self.stampp.mean
         #   y += Jy2muJy*self.stampp.median
         #   z = Jy2muJy*self.stampp.std
         #   print Jy2muJy*self.stampp.mean,Jy2muJy*self.stampp.median
            # Add it to the array - that's it!

    def writeSTD(self,running_std_S,m):
        """
        Write the running standard deviation to png
        """

        import numpy
        import matplotlib.pyplot as plt
        from math import sqrt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        xx=[m*el for el in range(len(running_std_S))]
        plt.plot(xx,running_std_S)
        x=numpy.arange(10,10**6)
        y=numpy.arange(10,10**6)
        y=[1/sqrt(el) for el in y]
        plt.plot(x,y)
        plt.xscale('log')
        plt.yscale('log')
        plt.title('std_S/muJy')
        plt.xlabel('sample')
        plt.ylabel('std_S/muJy')
        f='run_std_S.png'
        plt.savefig(f)
        print 'wrote %s' % f
        plt.close()

        return


    def writeFrame(self,frame,nn,f=None,flux=None):
        """
        Write latest mean frame to png
        """

        import numpy
        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        if flux=='S':
            cax=ax.imshow(frame,interpolation='nearest',vmin=-0.5,vmax=1.0)
        elif flux=='L':
            cax=ax.imshow(frame,interpolation='nearest')
        fig.colorbar(cax)
        plt.gca().invert_yaxis()
        if flux=='S':
            plt.title('mean_%s/muJy (std/muJy=%f, n=%i)' % (flux,numpy.std(frame),nn))
        elif flux=='L':
            plt.title('mean_%s/muJy (std/muJy=%e, n=%i)' % (flux,numpy.std(frame),nn))
        plt.xlabel('RA/pix')
        plt.ylabel('Dec/pix')
        if f is None: f='frames/frame_%s_%i.png' % (flux,nn)
        plt.savefig(f)
        print 'wrote %s' % f
        plt.close()

        return

    def pixels(self,radio_map,radio_header,noise_map,\
               f=None,simulate=False,report='report.txt',\
               noise_thresh=None,clip_thresh=None,directory=None,\
               upsampling_factor=None,masses='bon',\
               upsampled_map=None,upsampled_noise=None,\
               use_radio_posn=False,do_sim_posns=None,\
               scramble_posns=None,noise_only=None,\
               scramble_posns_seed=None,scramble_posns_mask=None):
        """
        For each K-band galaxy position,
        extract the radio pixel value at that position
        Return the fluxes, luminosities and stats (for use by a wrapper, e.g. wrapMzCut)
        noise_thresh in muJy and clip_thresh in sigma
        """
        import numpy
        from math import sqrt,log10
        from numpy.random import RandomState
        #from scipy.stats import morestats
        from astLib import astWCS
        import os
        from image_registration.fft_tools import upsample_image

        confidence_interval=0.68
        print 'Confidence interval is %4.2f' % confidence_interval

        if f is None: f='pixels.dat'
        #f='%s%i' % (f,self.seed)
        f='%s_%i.%s' % ('.'.join(f.split('.')[:-1]),self.seed,f.split('.')[-1])
        if directory is not None:
            f=os.path.join(directory,f)

        d=open(f,'w')
        #d.write('# n n_noise_clipped id alpha delta k x y z sovern noise weight S_muJy L_WperHz\n')
        d.write('# n n_noise_clipped id alpha delta k x y z sovern noise weight S_muJy L_WperHz M SFR_c SSFR_c SFR_y SSFR_y reliability mass_bon err_mass_bon sfr_bon err_sfr_bon sfr_gio err_sfr_gio mass_gio err_mass_gio\n')

        r=open(report,'a')

        if simulate: print 'simulating %i positions' % (len(self.gals))

        print 'Rejecting noise > %4.2f muJy and detections > %4.2e sigma' \
            % (noise_thresh,clip_thresh)

        n_ditched=0; n_detected=0
        self.Ss=numpy.array((1,1))
        self.Ns=numpy.array((1,1))
        self.Ws=numpy.array((1,1))
        self.Ls=numpy.array((1,1))
        self.zs=numpy.array((1,1))
        self.Ms=numpy.array((1,1))
        self.Ks=numpy.array((1,1))
        self.Kabs=numpy.array((1,1))      # For K-band absolute magnitude
        #self.Fs_c=numpy.array((1,1))      # For SFRs   Condon
        self.Rs_c=numpy.array((1,1))      # For SSFRs  Condon
        #self.logRs_c=numpy.array((1,1))   # For log SSFRs  Condon
        #self.Fs_y=numpy.array((1,1))      # For SFRs   Yun
        self.Rs_y=numpy.array((1,1))      # For SSFRs  Yun
        #self.logRs_y=numpy.array((1,1))   # For log SSFRs  Yun
        self.Fs_b=numpy.array((1,1))      # For SFRs Bonfield
        self.Rs_b=numpy.array((1,1))      # For SSFRs Bonfield
        self.Fs_g=numpy.array((1,1))      # For SFRs Giovannoli
        self.Rs_g=numpy.array((1,1))      # For SSFRs Giovannoli

        # For each galaxy:
        #WCS=astWCS.WCS('/Users/jtlz2/video/vla/VLA-VVDS-1.4GHZ.FITS')
        WCS=astWCS.WCS(radio_header,mode='pyfits')
        #test=open('test.txt','w')
        #test.write('# x y x1 y1\n')

        # Upsampling 16.4.14
        if upsampling_factor is not None:
            if upsampling_factor > 0:
                print '--> Beginning upsampling at x%i' % upsampling_factor
                upsampled_radio_map=upsample_image(radio_map[0,0,:,:],\
                                            upsample_factor=upsampling_factor,\
                                            use_numpy_fft=False)
                print '--> Finished upsampling radio map'
                print '--> Beginning upsampling at x%i' % upsampling_factor
                upsampled_noise_map=upsample_image(noise_map[0,0,:,:],\
                                            upsample_factor=upsampling_factor,\
                                            use_numpy_fft=False)
                print '--> Finished upsampling noise map'
            elif upsampling_factor < 0:
                print '--> Reading upsampled radio map from file %s'%upsampled_map
                radioArray=pyfits.open(upsampled_map)
                upsampled_radio_map=radioArray[0].data
                print '--> Reading upsampled noise map from file %s'%upsampled_noise
                noiseArray=pyfits.open(upsampled_noise)
                upsampled_noise_map=noiseArray[0].data
                if upsampling_factor is not None:
                    upsampling_factor=abs(upsampling_factor)
        else:
            pass

        # Read simulated positions from file, if required
        if do_sim_posns is not None:
            posns=numpy.loadtxt(os.path.join(do_sim_posns))
            nposns=len(posns[:,0])
            ngals=len(self.gals)
            assert(nposns<=ngals), 'Too few hosts! %i v. %i' % (nposns,ngals)
            self.gals=self.gals[:nposns]

        if noise_only is not None:
            print 'Drawing fluxes from a gaussian of width %f uJy' % noise_only
            numpy.random.seed(seed=9794)
            seeds=numpy.random.normal(0.0,noise_only,len(self.gals))

        if scramble_posns:
            numpy.random.seed(seed=scramble_posns_seed)
            # Hack:
            # Attempt to mask out positions with known sources
            if scramble_posns_mask:
                #dat=numpy.genfromtxt('/Users/jtlz2/video/20131007/VIDEO-cfhtls-D1_2013-10-07_fullcat_spec_fluxes_w_radio_jz_for_bayestack_140526.dat')
                #dat=dat[numpy.where(dat[:,37].astype(numpy.float64)>-99.0)]
                maskf='/Users/jtlz2/video/vla/Bondi2003_posns.txt'
                print 'Fetching mask positions from %s' % maskf
                dat=numpy.genfromtxt(maskf)
                #print dat
                # Identify positions to mask
                #xys=[astWCS.WCS.wcs2pix(WCS,gal.alpha,gal.delta) for gal in self.gals if gal.reliability >=0]
                #xys=[astWCS.WCS.wcs2pix(WCS,dat[iS,38],dat[iS,39]) for iS in range(len(dat))]
                xys=[astWCS.WCS.wcs2pix(WCS,dat[iS,0],dat[iS,1]) for iS in range(len(dat))]

                print len(xys)
                #numpy.savetxt('xys.txt',xys)
                xxyys=[(int(numpy.round(x*float(upsampling_factor))),int(numpy.round(y*float(upsampling_factor)))) for (x,y) in xys]
                coords = zip(*xxyys)
                #print upsampled_radio_map.shape,len(coords),len(coords[0])
                numpy.savetxt('xxyys.txt',xxyys)
                mask = sparse.coo_matrix((numpy.ones(len(coords[0])),coords),shape=upsampled_radio_map.shape,dtype=bool)
                #print mask
                mask2=mask.toarray()
                print numpy.count_nonzero(mask2)
                #mask2=dilation(mask3,disk(4*upsampling_factor)).astype(bool)
                #struct = scipy.ndimage.generate_binary_structure(mask3.ndim,10*upsampling_factor)
                #struct=tophat(1.0,21*upsampling_factor,10*upsampling_factor)[1]
                #mask2=scipy.ndimage.morphology.binary_dilation(mask3,structure=struct).astype(bool)
                #gD=upsampled_radio_map.shape[0]
                SURVEY_NOISE=16.2*Jy2muJy # hack
                radioSynthBeamFWHM=4.0 # pix
                maskRadius=5.0 # was 5.0
                rad=maskRadius*radioSynthBeamFWHM*upsampling_factor
                gD2=maskRadius*radioSynthBeamFWHM*upsampling_factor
                struct=tophat(1.0,(2*gD2)+1,rad)[1]
                frame=50*upsampling_factor
                accumulator=numpy.zeros(upsampled_radio_map.shape)
                #numpy.pad(accumulator,gD2)
                for iS in range(len(xxyys)):
                    #print iS,mask.row[iS],mask.col[iS]
                    #y,x = numpy.ogrid[-mask.row[iS]:gD-mask.row[iS],-mask.col[iS]:gD-mask.col[iS]]
                    #mask2=numpy.bitwise_or(mask2,x*x + y*y <= r*r)
                    #print accumulator[mask.col[iS]-gD2:mask.col[iS]+gD2+1,mask.row[iS]-gD2:mask.row[iS]+gD2+1]
                    #print struct
                    #xc=mask.row[iS]; yc=mask.col[iS]
                    #xc=10; yc=10
                    #ytop,xtop=upsampled_radio_map.shape
                    #xmin=max(xc-rad,0); xmax=min(xc+rad+1,xtop)
                    #ymin=max(yc-rad,0); ymax=min(yc+rad+1,ytop)
                    print iS+1#,xmin,xmax,ymin,ymax,xc,yc
                    accumulator[mask.col[iS]-gD2:mask.col[iS]+gD2+1,mask.row[iS]-gD2:mask.row[iS]+gD2+1]+=struct
                    #print accumulator[ymin:ymax,xmin:xmax].shape,struct.shape
                    #dya,dxa=accumulator[ymin:ymax,xmin:xmax].shape
                    #print struct[:200,:200]
                    #accumulator[ymin:ymax,xmin:xmax]+=struct[gD2-dya//2:gD2+1+dya//2,gD2-dxa//2:gD2+1+dxa//2]
                    #sys.exit(0)
                mask3=numpy.ma.make_mask(accumulator)
                #print mask3,type(mask3)
            ##    mask4=numpy.where(upsampled_radio_map>5.0*SURVEY_NOISE)
                #pyfits.writeto('accum.fits',accumulator,clobber=True)
                #pyfits.writeto('mask.fits',mask3.astype(numpy.int64),clobber=True)
                del accumulator
                print mask3.size,numpy.count_nonzero(mask3)
                print mask3[frame:upsampled_radio_map.shape[0]-1-frame,frame:upsampled_radio_map.shape[1]-1-frame].size,numpy.count_nonzero(mask3[frame:upsampled_radio_map.shape[0]-1-frame,frame:upsampled_radio_map.shape[1]-1-frame])
                assert(len(self.gals)<numpy.sum(~mask3)), \
                  '%i draws > %i free pix' % (len(self.gals),numpy.sum(~mask3))

                #import itertools
                #sparse_fluxes=Jy2muJy*numpy.random.choice(upsampled_radio_map[~mask.toarray()].flatten(),size=len(self.gals),replace=False)
                #indices=list(numpy.ndindex(upsampled_radio_map.shape))
                #nzm=numpy.nonzero(mask3)
                #indices_to_mask=itertools.izip(nzm[0].ravel(),nzm[1].ravel())
                #print len(indices),len(indices_to_mask)
                #indices_remaining=(i for i in xrange(numpy.ndindex(upsampled_radio_map.shape)) if i not in indices_to_mask)
                #print len(indices_remaining)
                #draws=numpy.random.randint(0,len(indices_remaining),len(self.gals))
                #ir=(indices_remaining[i] for i in draws)
                #sparse_fluxes=Jy2muJy*numpy.random.choice(upsampled_radio_map[~mask3].flatten(),size=len(self.gals),replace=False)
                #sparse_flux_indices=numpy.random.choice(numpy.arange())
                ndraws=len(self.gals)
                ni=0; coords=numpy.zeros((ndraws,2))
                sparse_fluxes=numpy.zeros(ndraws)
                while (ni<ndraws):
                    x=numpy.random.randint(frame,upsampled_radio_map.shape[0]-1-frame)
                    y=numpy.random.randint(frame,upsampled_radio_map.shape[1]-1-frame)
                    if not mask3[y,x]:
                        coords[ni,:]=[y,x]
                        #print '%12i   x: %5i y: %5i' % (large_table[x][y], x,y)
                        sparse_fluxes[ni]=Jy2muJy*upsampled_radio_map[y,x]
                        #if sparse_fluxes[ni]>85.0:
                        #    print ni,sparse_fluxes[ni],x,y,mask3[y-10:y+10,x-10:x+10]
                        #    print upsampled_radio_map[y-10:y+10,x-10:x+10]
                        #    sys.exit(0)
                        mask3[y,x]=True
                        ni+=1
                pyfits.writeto('mask.fits',mask3.astype(numpy.int64),clobber=True)
                numpy.savetxt('fluxes-masked.txt',sparse_fluxes)
                del mask2,mask3
                print 'Look in %s' % 'fluxes-masked.txt'
                #sys.exit(0)

        for n,gal in enumerate(self.gals):
            #if gal.reliability>=0: print gal.reliability
            if simulate or scramble_posns:
                if scramble_posns:
                    x=numpy.random.randint(0,radio_header['NAXIS1'])
                    y=numpy.random.randint(0,radio_header['NAXIS2'])
                else:
                    x=self.prng.randint(0,radio_header['NAXIS1'])
                    y=self.prng.randint(0,radio_header['NAXIS2'])
                xx=x; yy=y
            else:
                # wcs2pix returns floats numbered from 0.0,0.0 [0,0]
                if not use_radio_posn:
                    x,y=astWCS.WCS.wcs2pix(WCS,gal.alpha,gal.delta)
                else:
                    x,y=astWCS.WCS.wcs2pix(WCS,gal.alpha_radio,gal.delta_radio)

            if do_sim_posns:
                x=posns[n,1]-1; y=posns[n,2]-1 # Column 0 is just an index

            # Now round these to nearest integer
            xx=x; yy=y # 15.4.14 Make a copy
            x=int(numpy.round(x))
            y=int(numpy.round(y))
            xxx=xx; yyy=yy

                #line='%i %i %i %i\n'%(x,y,x1,y1)
                #test.write(line)

            # Exclude regions on the radio image with sigma_1.4 > 40 muJy
            #noise=Jy2muJy*noise_map[0,0,y-1,x-1]
            if upsampling_factor is not None:
                xx=int(numpy.round(xx*float(upsampling_factor)))
                yy=int(numpy.round(yy*float(upsampling_factor)))
                noise=Jy2muJy*upsampled_noise_map[yy,xx]
            else:
                noise=Jy2muJy*noise_map[0,0,y,x]

            weight=1.0/(noise*noise)
            #print noise,noise_thresh
            if (noise > noise_thresh):
                n_ditched+=1
                continue

            # Extract the pixel fluxes and calculate the luminosities
            # Upsampling 16.4.14
            if upsampling_factor is not None:
                if not scramble_posns_mask:
                    S=Jy2muJy*upsampled_radio_map[yy,xx] # Conversion to int already done ^^^
                else:
                    S=sparse_fluxes[n]
                    #S=Jy2muJy*upsampled_radio_map[ir[n]]
                    yy,xx=coords[n,:]
        ##      S2=Jy2muJy*radio_map[0,0,y,x]
                #print xx,yy,S
                #sys.exit(0)
                fwhm=4.0*abs(upsampling_factor) #3.8148
                sig=fwhm/2.35482 # FWHM = 6" x 6" = 4 x 4 pixels
                #sig=1.69
                amp=1.0/(2.0*pi*sig*sig)
                PSF=GaussianPSF(float(sig),amplitude=1.0)
                fac=pi/(4.0*log(2.0))

                #S4 = Jy2muJy*psf_photometry(upsampled_radio_map[:,:],[(xx,yy)],PSF)[0]/fac/fwhm/fwhm
                #Starget=1.0e3*gal.radio_pflux
                #S4=0.0; dScurrent=0.0; dSclosest=1.0e10; xS=0; yS=0
                ##a[numpy.sum(numpy.square(numpy.abs(a-a0)),1).argmin()]
                #for ix in xrange(-2*abs(upsampling_factor),3*abs(upsampling_factor)):
                #    for iy in xrange(-2*abs(upsampling_factor),3*abs(upsampling_factor)):
                #        Swalk=Jy2muJy*psf_photometry(upsampled_radio_map[:,:],[(xx+0.5+ix,yy+0.5+iy)],PSF)[0]/fac/fwhm/fwhm
                #        dScurrent=abs(Swalk-Starget)
                #        if dScurrent<dSclosest:
                #            dSclosest=dScurrent
                #            S4=Swalk
                #            xS=ix; yS=iy
                #        print gal.id,xx,yy,ix,iy,dScurrent,dSclosest,Starget,S4,Swalk
                #print gal.id,xS,yS,S,Starget,S4,S4/1.0e3-gal.radio_pflux,100.0*(S4/1.0e3-gal.radio_pflux)/gal.radio_pflux,'\n'
                #print gal.id,S,S4,S4-S,100.0*(S4-S)/S
                #print n,xxx,yyy,x,y,xx/float(upsampling_factor),\
                #  yy/float(upsampling_factor),xx,yy,\
                #  x*upsampling_factor,y*upsampling_factor,\
                #  S,S2,100.0*(S-S2)/S2,S4,100.0*(S-S4)/S4
                # Aperture photometry not applied/implemented for upsampled calculations
                #noise4=Jy2muJy*psf_photometry(upsampled_noise_map[:,:],[(xx,yy)],PSF)[0]/fac/fwhm/fwhm
                #norm=psf_photometry(numpy.ones((int(2.0*fwhm)+1,int(2.0*fwhm)+1)),[(7,7)],PSF)[0]/fac/fwhm/fwhm
                #noise4 *= 1.0/norm
                # S=S4
            else:
                S=Jy2muJy*radio_map[0,0,y,x] # ***** NB Arrays count from zero!!
                fwhm=4.0 #3.8148
                sig=fwhm/2.35482 # FWHM = 6" x 6" = 4 x 4 pixels
                #sig=1.69
                amp=1.0/(2.0*pi*sig*sig)
                PSF=GaussianPSF(float(sig),amplitude=1.0)
                fac=pi/(4.0*log(2.0))
                #print x,xx,y,yy
                #S4 = Jy2muJy*psf_photometry(radio_map[0,0,:,:],[(x+0.5,y+0.5)],PSF)[0]/fac/fwhm/fwhm
                #print x,y,S,S4
                #raw_input('')
                #print numpy.shape(radio_map[0,0,yy-6:yy+6,xx-6:xx+6])
                #S5 = Jy2muJy*sum2DGaussian(radio_map[0,0,y-7:y+7,x-7:x+7],sig,1.0,yy,xx)
                                #S4 *= 3.269
                #noise4=Jy2muJy*psf_photometry(noise_map[0,0,:,:],[(xx,yy)],PSF)[0]/fac/fwhm/fwhm

                #Starget=1.0e3*gal.radio_pflux
                #S4=0.0; dScurrent=0.0; dSclosest=1.0e10
                #for ix in xrange(-2,3):
                #    for iy in xrange(-2,3):
                #        Swalk=Jy2muJy*psf_photometry(radio_map[0,0,:,:],[(x+0.5+ix,y+0.5+iy)],PSF)[0]/fac/fwhm/fwhm
                #        dScurrent=abs(Swalk-Starget)
                #        if dScurrent<dSclosest:
                #            dSclosest=dScurrent
                #            S4=Swalk
                #        #print ix,iy,dScurrent,dSclosest,S4,Swalk

                #        noise4=Jy2muJy*psf_photometry(noise_map[0,0,:,:],[(x+0.5,y+0.5)],PSF)[0]/fac/fwhm/fwhm
                #noise5=Jy2muJy*numpy.mean(noise_map[0,0,xx-fwhm:xx+fwhm,yy-fwhm:yy+fwhm])
                # Normalize by the volume of the gaussian kernel
                lside=int(8.0*sig)+1; icentre=int(lside/2.0)
                #norm=psf_photometry(numpy.ones((lside,lside)),[(icentre,icentre)],PSF)[0]/fac/fwhm/fwhm
                #print lside,icentre,norm
                #noise4 *= 1.0/norm
                #print gal.id,S,Starget,S4,S4/1.0e3-gal.radio_pflux,100.0*(S4/1.0e3-gal.radio_pflux)/gal.radio_pflux
                #sys.exit(0)
                #noise=1.0/sqrt(weight)
                #noise4=Jy2muJy*sqrt(weight4)/fac/fwhm/fwhm#/14.0
                #print gal.id,noise,noise4,noise4-noise,noise4/noise,100.0*(noise4-noise)/noise
                #print weight,weight4,weight4/weight,100.0*(weight4-weight)/weight
                #print
                #noise,noise4,noise4-noise,noise4/noise,100.0*(noise4-noise)/noise#,S5,S5-S,100.0*(S5-S)/S,S5,S5-S4,100.0*(S5-S4)/S4
                
                #S=S4; noise=noise4

            # Bilinear interpolation 15.4.14
            #x1=int(xx); x2=int(xx+1); y1=int(yy); y2=int(yy+1)
            #x1=numpy.floor(xx); x2=numpy.ceil(xx); y1=numpy.floor(yy); y2=numpy.ceil(yy)
            #p1=(x1,y1,radio_map[0,0,y1,x1])
            #p2=(x1,y2,radio_map[0,0,y2,x1])
            #p3=(x2,y1,radio_map[0,0,y1,x2])
            #p4=(x2,y2,radio_map[0,0,y2,x2])
            #S=Jy2muJy*bilinear_interpolation(xx,yy,[p1,p2,p3,p4])

            # Optionally just calculate the flux to be a gaussian draw
            #     from SURVEY_NOISE
            if noise_only is not None:
                #print 'Drawing fluxes from SURVEY_NOISE'
                S=seeds[n]
            
            # Exclude pixels with |S/N| >= 5 (Dunne) ***ONLY FOR MEAN***
            # It would be nice to have allowed this in // with the median,
            #    but probably too late now
            #sovern=radio_map[0,0,y-1,x-1]/noise_map[0,0,y-1,x-1]
            #sovern=radio_map[0,0,y,x]/noise_map[0,0,y,x]
            sovern=S/noise
            if sovern >= clip_thresh or sovern <= -clip_thresh:
                n_detected+=1
                continue

            self.Ns=numpy.append(self.Ns,noise)
            self.Ws=numpy.append(self.Ws,weight)

            #pyfits.writeto('upsampled.fits',upsampled_map)
            self.Ss=numpy.append(self.Ss,S)
            L=gal.factor*S/Jy2muJy
            self.Ls=numpy.append(self.Ls,L)

            # Collate redshifts, stellar masses and K mags
            z=gal.z_best
            self.zs=numpy.append(self.zs,z)
            if masses=='bon':
                M=gal.mass_bon
            elif masses=='gio':
                M=gal.mass_gio     
            self.Ms=numpy.append(self.Ms,M)
            K=gal.kmag
            self.Ks=numpy.append(self.Ks,K)
            Kabs=gal.kabsmag
            self.Kabs=numpy.append(self.Kabs,Kabs)

            # Calculate Condon SFR and SSFR
            SFR_c=1.2006e-21*L
            #self.Fs_c=numpy.append(self.Fs_c,SFR_c)
            SSFR_c=SFR_c/(10**M) # Since M is in units of log_10(M)
            self.Rs_c=numpy.append(self.Rs_c,SSFR_c)
            #log_SSFR_c=log10(SFR_c)-M
            #self.logRs_c=numpy.append(self.logRs_c,log_SSFR_c)

            # Calculate Yun SFR and SSFR
            SFR_y=5.9e-22*L
            #self.Fs_y=numpy.append(self.Fs_y,SFR_y)
            SSFR_y=SFR_y/(10**M) # Since M is in units of log_10(M)
            self.Rs_y=numpy.append(self.Rs_y,SSFR_y)
            #log_SSFR_y=log10(SFR_y)-M
            #self.logRs_y=numpy.append(self.logRs_y,log_SSFR_y)

            # Collate Bonfield SFR and SSFR
            SFR_bon=10**gal.sfr_bon
            self.Fs_b=numpy.append(self.Fs_b,SFR_bon)
            SSFR_bon=SFR_bon/(10**gal.mass_bon) # Since M is in units of log_10(M)
            self.Rs_b=numpy.append(self.Rs_b,SSFR_bon)

            # Collate Giovannoli SFR and SSFR
            SFR_gio=10**gal.sfr_gio
            self.Fs_g=numpy.append(self.Fs_g,SFR_gio)
            SSFR_gio=SFR_gio/(10**gal.mass_gio) # Since M is in units of log_10(M)
            self.Rs_g=numpy.append(self.Rs_g,SSFR_gio)


            #w='%6i %6i %10i %f %f %f %4i %4i %4.2f %4.2f %4.2f %e %4.2f %e' \
            #    % (n,(n_ditched+n_detected),gal.id,gal.alpha,gal.delta,gal.kmag,\
            #       x,y,gal.z_best,sovern,noise,weight,S,L)
            w='%6i %6i %10i %12.9f %12.9f %f %4f %4f %5.3f %5.3f %5.3f %e %5.3f %e %f %e %e %e %e %f %f %f %e %e %f %f %e %e' \
                % (n,(n_ditched+n_detected),gal.id,gal.alpha,gal.delta,gal.kmag,\
                   xx,yy,gal.z_best,sovern,noise,weight,S,L,M,SFR_c,SSFR_c,SFR_y,\
                   SSFR_y,gal.flag,gal.mass_bon,gal.err_mass_bon,gal.sfr_bon,\
                   gal.err_sfr_bon,gal.sfr_gio,gal.err_sfr_gio,gal.mass_gio,\
                   gal.err_mass_gio)

            d.write('%s\n'%w)
            if n > 0 and n % 10000 == 0:
                if not simulate: print w
                d.flush()

        d.close()


        # Now that we have per-galaxy values, generate summary stats
        # Initialize stats dict
        stats={}

        # The first two values are 1.0 - discard.
        # I did know why this was but have forgotten..
        Ss=self.Ss[2:]; Ls=self.Ls[2:]; zs=self.zs[2:]; Kabs=self.Kabs[2:]
        Ns=self.Ns[2:]; Ws=self.Ws[2:]; Ms=self.Ms[2:]; Ks=self.Ks[2:]
        #Fs_c=self.Fs_c[2:]; Fs_y=self.Fs_y[2:]
        Rs_c=self.Rs_c[2:]; Rs_y=self.Rs_y[2:]
        #logRs_c=self.Rs_c[2:]; logRs_y=self.Rs_y[2:]
        Fs_b=self.Fs_b[2:]; Rs_b=self.Rs_b[2:]
        Fs_g=self.Fs_g[2:]; Rs_g=self.Rs_g[2:]

        if numpy.size(Ss) <= 2:
            stats['nsrc']=0
            return Ss,Ls,Ws,Ms,Ks,Kabs,stats

        # Calculate error bars on S's
        # The confidence interval *is* centred on the median
        # [0] is for mean rather than var or std
        S_err_med_low,S_err_med_high,S_error_median=err(Ss,alpha=confidence_interval)

        # Start the report to screen and to file f
        line= '==\n== Report for %s:\n==' % f
        print line; r.write('%s\n'%line)

        # Report flux statistics S -- mean and median
        line= \
            '== S / muJy  mean %5.3f       +/- %5.3f [std %5.3f] median <%5.3f> +/- %5.3f'\
            % (numpy.average(Ss), numpy.std(Ss)/sqrt(len(Ss)-1),\
            numpy.std(Ss), numpy.median(Ss), S_error_median)
        print line; r.write('%s\n'%line)
        stats['S_median']     = numpy.median(Ss)
        stats['S_median_err'] = S_error_median
        stats['S_median_err_low'] = S_err_med_low
        stats['S_median_err_high'] = S_err_med_high

        # Report flux statistics S -- weighted mean
        (dummy,weighted_std_S,weights_sum_S)=weighted_avg_and_std(Ss,Ws)
        line= \
            '== S / muJy wmean %5.3f       +/- %5.3f [std %5.3f] {sum_w %5.3f}'\
            % (numpy.average(Ss,weights=Ws),weighted_std_S/sqrt(len(Ss)-1),\
               weighted_std_S,weights_sum_S)
        print line; r.write('%s\n'%line)

        # Calculate luminosity statistics L -- mean and median -- but not for sims
        if not simulate:

            # Calculate median L and report
            stats['L_median']     = numpy.median(Ls)
            stats['L_median_err_low'],stats['L_median_err_high'],stats['L_median_err'] \
                = err(Ls,alpha=confidence_interval)

            line= '== L / W/Hz  mean %5.3e   +/- %5.3e   median <%5.3e> +/- %5.3e' \
                % (numpy.average(Ls), numpy.std(Ls)/sqrt(len(Ls)-1),\
                    stats['L_median'],stats['L_median_err'])
            print line; r.write('%s\n'%line)

            # Report luminosity statistics L -- weighted mean
            (dummy,weighted_std_L,weights_sum_L)=weighted_avg_and_std(Ls,Ws)
            line= \
                '== L / W/Hz wmean %5.3e   +/- %5.3e   [std %5.3e] {sum_w %5.3f}'\
                % (numpy.average(Ls,weights=Ws),weighted_std_L/sqrt(len(Ls)-1),\
                   weighted_std_L,weights_sum_L)
            print line; r.write('%s\n'%line)

            # Calculate median Condon SFR and report
            stats['SFR_condon_median']     = 1.2006e-21*stats['L_median']
            stats['SFR_condon_median_err'] = 1.2006e-21*stats['L_median_err']
            stats['SFR_condon_median_err_low'] = 1.2006e-21*stats['L_median_err_low']
            stats['SFR_condon_median_err_high'] = 1.2006e-21*stats['L_median_err_high']

            line= '==  SFR (Condon; median)    / M_sun/yr <%5.3f> +/- %5.3f' % \
                (stats['SFR_condon_median'], stats['SFR_condon_median_err'])
            print line; r.write('%s\n'%line)

            # Calculate median Yun SFR and report
            SFR_yun_median=5.9e-22*stats['L_median']
            # Work out the error for SFR_yun, using error given for the conversion
            # Add these in quadrature
            A=(1.8e-22/5.9e-22)
            B=(stats['L_median_err']/stats['L_median'])
            SFR_yun_error_median=SFR_yun_median*sqrt(A**2+B**2)
            stats['SFR_yun_median']     = SFR_yun_median
            stats['SFR_yun_median_err'] = SFR_yun_error_median
            line= '==  SFR (Yun;    median)    / M_sun/yr <%5.3f> +/- %5.3f' % \
                (stats['SFR_yun_median'], stats['SFR_yun_median_err'])
            print line; r.write('%s\n'%line)

            # Calculate median M and report
            stats['M_median']     = numpy.median(Ms)
            stats['M_median_err_low'],stats['M_median_err_high'],stats['M_median_err'] \
                = err(Ms,alpha=confidence_interval)
            
            line= '== M_st                     / M_sun <%5.3f>     +/- %5.3f' %\
                (stats['M_median'],stats['M_median_err'])
            print line; r.write('%s\n'%line)

            # Calculate median Condon SSFR and report
            stats['SSFR_condon_median']     = 1.0e9*numpy.median(Rs_c) # /yr -> /Gyr
            stats['SSFR_condon_median_err_low'],stats['SSFR_condon_median_err_high'],\
                stats['SSFR_condon_median_err'] \
                = err(1.0e9*Rs_c,alpha=confidence_interval)

            line= '== SSFR (Condon; median)    / M_sun/yr <%5.3f> +/- %5.3f' % \
                (stats['SSFR_condon_median'],stats['SSFR_condon_median_err'])
            print line; r.write('%s\n'%line)

            # Calculate median Condon log_SSFR and report
            # THIS FAILS because sometimes SFR_c (Fs_c) is negative
            # INSTEAD one needs to take the other route
            #logRs_c=numpy.log10(Fs_c)-Ms
            if stats['SSFR_condon_median'] > 0.0:
                stats['log_SSFR_condon_median'] = log10(stats['SSFR_condon_median'])
            else:
                stats['log_SSFR_condon_median'] = float('nan')
            stats['log_SSFR_condon_median_err'] = \
                leuler*stats['SSFR_condon_median_err']/stats['SSFR_condon_median']
            line= '== logSSFR (Condon; median) / M_sun/yr <%5.3f> +/- %5.3f' % \
                (stats['log_SSFR_condon_median'],stats['log_SSFR_condon_median_err'])
            print line; r.write('%s\n'%line)

            # Calculate median Yun SSFR and report
            stats['SSFR_yun_median']     = 1.0e9*numpy.median(Rs_y) # /yr -> /Gyr
            stats['SSFR_yun_median_err'] = err(1.0e9*Rs_y,alpha=confidence_interval)[2]
            line= '== SSFR (Yun;    median)    / M_sun/yr <%5.3f> +/- %5.3f [uninflated]' % \
                (stats['SSFR_yun_median'],stats['SSFR_yun_median_err'])
            print line; r.write('%s\n'%line)

            # Calculate median Yun log_SSFR and report
            # THIS FAILS because sometimes SFR_y (Fs_y) is negative
            # INSTEAD one needs to take the other route
            #logRs_y=numpy.log10(Fs_y)-Ms
            if stats['SSFR_yun_median'] > 0.0:
                stats['log_SSFR_yun_median'] = log10(stats['SSFR_yun_median'])
            else:
                stats['log_SSFR_yun_median'] = float('nan')
            stats['log_SSFR_yun_median_err'] = \
                leuler*stats['SSFR_yun_median_err']/stats['SSFR_yun_median']

            line= '== logSSFR (Yun;    median) / M_sun/yr <%5.3f> +/- %5.3f' % \
                (stats['log_SSFR_yun_median'],stats['log_SSFR_yun_median_err'])
            print line; r.write('%s\n'%line)


            # Calculate median Bonfield SSFR and SFR and report
            stats['SFR_bon_median']      = numpy.median(Fs_b) # /yr
            stats['SFR_bon_median_err']  = err(Fs_b,alpha=confidence_interval)[2]
            stats['SSFR_bon_median']     = 1.0e9*numpy.median(Rs_b) # /Gyr -?-> /yr
            stats['SSFR_bon_median_err'] = err(1.0e9*Rs_b,alpha=confidence_interval)[2]
            if stats['SSFR_bon_median'] > 0.0:
                stats['log_SSFR_bon_median'] = log10(stats['SSFR_bon_median'])
            else:
                stats['log_SSFR_bon_median'] = float('nan')
            stats['log_SSFR_bon_median_err'] = \
                leuler*stats['SSFR_bon_median_err']/stats['SSFR_bon_median']

            # WORKING HERE
            # Calculate median Giovannoli SSFR and SFR and report
            stats['SFR_gio_median']      = numpy.median(Fs_g) # /yr
            stats['SFR_gio_median_err']  = err(Fs_g,alpha=confidence_interval)[2]
            stats['SSFR_gio_median']     = 1.0e9*numpy.median(Rs_g) # /Gyr -?-> /yr
            stats['SSFR_gio_median_err'] = err(1.0e9*Rs_g,alpha=confidence_interval)[2]
            if stats['SSFR_gio_median'] > 0.0:
                stats['log_SSFR_gio_median'] = log10(stats['SSFR_gio_median'])
            else:
                stats['log_SSFR_gio_median'] = float('nan')
            stats['log_SSFR_gio_median_err'] = \
                leuler*stats['SSFR_gio_median_err']/stats['SSFR_gio_median']

            if masses=='bon':
                line= '== logSSFR (Bon;    median) / M_sun/yr <%5.3f> +/- %5.3f' % \
                (stats['log_SSFR_bon_median'],stats['log_SSFR_bon_median_err'])
            elif masses=='gio':
                line= '== logSSFR (Gio;    median) / M_sun/yr <%5.3f> +/- %5.3f' % \
                (stats['log_SSFR_gio_median'],stats['log_SSFR_gio_median_err'])
            print line; r.write('%s\n'%line)

            # Calculate median K...
            stats['K_median']     = numpy.median(Ks)
            stats['K_median_err'] = err(Ks,alpha=confidence_interval)[2]

            # Calculate median Kabs...
            stats['Kabs_median']     = numpy.median(Kabs)
            stats['Kabs_median_err'] = err(Kabs,alpha=confidence_interval)[2]

            # Calculate median z...
            stats['z_median']     = numpy.median(zs)
            stats['z_median_err'] = err(zs,alpha=confidence_interval)[2]

            # ... and report median z, M_stellar, K and Kabs
            line= '== median z <%5.3f>  M_st <%5.3f>  K <%5.3f>  Kabs <%5.3f>' \
                % (stats['z_median'],stats['M_median'],stats['K_median'],stats['Kabs_median'])
            print line; r.write('%s\n'%line)

        # Calculate fraction detected and noisy and report
        frac_detected=float(n_detected)/float(n)
        frac_ditched=float(n_ditched)/float(n)
        stats['frac_detected']=frac_detected
        stats['frac_ditched']=frac_ditched
        stats['nsrc']=n+1
        line= '== numbers: %i in bin, o/w detected %5.3f noisy %5.3f per cent \n==' \
            % (n+1,100.0*frac_detected,100.0*frac_ditched)
        print line; r.write('%s\n'%line)

        r.close()

        return Ss,Ls,Ws,Ms,Ks,Kabs,stats

#-----------------------------------------------------------

class squash(object):
    """
    This is a copy of the final postage stamps,
    to avoid having to regenerate each of many cutouts
    every time the science methods change
    """

    def __init__(self,cube):
        #self.cube=cube
        self.stamps=cube.stamps
        self.stamps_L=cube.stamps_L
        self.gals=cube.gals
        self.nx=cube.nx
        self.ny=cube.ny
        self.simulate=cube.simulate
        self.collapsed=False
        self.corner_noise=False
        self.shuffle=False

    def collapse(self,average=None,corner_noise=False):
        """
        This extracts science results (for mean or median) from the final postage stamp
        A bit slow
        """

        import numpy
        import matplotlib.pyplot as plt
        import copy
        from math import sqrt

        if average==None:
            print 'Choice of average required (mean or median)!'
            return

        self.corner_noise=corner_noise

        print 'Cube has dimensions', numpy.shape(self.stamps)

        #print stamp.x,stamp.y,stamp.gal.alpha,stamp.gal.delta

        # Stack the cutouts in flux
        self.average=average
        if average=='mean':
            self.stack_avg=numpy.average(self.stamps,axis=0)
        elif average=='median':
            self.stack_avg=numpy.median(self.stamps,axis=0)

        if self.corner_noise:
            self.stack_noise=numpy.std(self.stack_avg[:5,:5])
        else:
            self.stack_noise=numpy.std(self.stack_avg[:,:])


        # Report S_peak
        if self.simulate:
            # Assess simulated maps:
            if self.nx==21:
                x_max,y_max=(9,9) # <------- watch out for this later!
            elif self.ny==41:
                x_max,y_max=(19,19) # <------- watch out for this later!
            print '(For simulation, taking peak as x=y=%i)' % x_max
            print 'Peak S %6.3f (+/- %6.3f) muJy at %i,%i' % \
                (Jy2muJy*self.stack_avg[x_max,y_max],Jy2muJy*self.stack_noise,x_max,y_max)
        else:
            x_max,y_max=numpy.unravel_index(self.stack_avg.argmax(), self.stack_avg.shape)
            print 'Peak S %6.3f (+/- %6.3f) muJy at %i,%i' % \
                (Jy2muJy*numpy.max(self.stack_avg),Jy2muJy*self.stack_noise,x_max,y_max)

        self.x_max_S=x_max
        self.y_max_S=y_max

        # Sanity check on position of peak
        if (x_max > self.nx/2+4 or x_max < self.nx/2-4\
            or y_max > self.ny/2+4 or y_max < self.ny/2-4):
            print '**Warning: maximum is not near centre of image'

        # Write the peak fluxes (for each set of cutouts) to file,
        # along the line of cutouts
        # These can be checked for their mean and median values
        if self.simulate:
            f=open('avg_sim.txt','w')
        else:
            f=open('avg_dat.txt','w')
        for el in range(len(self.stamps[:,self.x_max_S,self.y_max_S])):
            f.write('%f\n' % self.stamps[el,self.x_max_S,self.y_max_S])
        f.close()

        # Stack the cutouts in luminosity
        if average=='mean':
            self.stack_avg_L=numpy.average(self.stamps_L,axis=0)
        elif average=='median':
            self.stack_avg_L=numpy.median(self.stamps_L,axis=0)

        if self.corner_noise:
            self.stack_noise_L=numpy.std(self.stack_avg_L[:5,:5])
        else:
            self.stack_noise_L=numpy.std(self.stack_avg_L[:,:])

        # Report L
        if self.simulate:
            if self.nx==21:
                x_max,y_max=(9,9)
            elif self.ny==41:
                x_max,y_max=(19,19) # <------- watch out for this later!
            print '(For simulation, taking peak as x=y=%i)' % x_max
            print 'Peak L %6.2e (+/- %6.2e) W/Hz at %i,%i' % (self.stack_avg_L[x_max,y_max],\
                                                              self.stack_noise_L,\
                                                              x_max,y_max)
        else:
            x_max,y_max=numpy.unravel_index(self.stack_avg_L.argmax(),\
                                            self.stack_avg_L.shape)
            print 'Peak L %6.2e (+/- %6.2e) W/Hz at %i,%i' % (numpy.max(self.stack_avg_L),\
                                                              self.stack_noise_L,\
                                                              x_max,y_max)
        self.x_max_L=x_max
        self.y_max_L=y_max

        # Sanity check on position of peak
        if (x_max > self.nx/2+4 or x_max < self.nx/2-4\
            or y_max > self.ny/2+4 or y_max < self.ny/2-4):
            print '**Warning: maximum is not near centre of image'

        # Report SFRs
        print '--> SFR_dunne M_solar/year %7.3f +/- %7.3f' % \
            (1.2006e-21*numpy.max(self.stack_avg_L),\
            1.2006e-21*self.stack_noise_L)

        # Work out the error for SFR_yun
        # **** THIS MUST BE CHECKED CAREFULLY
        print '**** WARNING: !!The next line must be checked carefully!!'
        A=(1.8e-22/5.9e-22)
        B=(self.stack_noise_L/numpy.max(self.stack_avg_L))
        quad_error_yun=(5.9e-22*numpy.max(self.stack_avg_L))*sqrt(A**2+B**2)
        print '--> SFR_yun M_solar/year %7.3f +/- %7.3f' % \
            (5.9e-22*numpy.max(self.stack_avg_L), quad_error_yun)

        self.collapsed=True

    def noiseHistogram(self):
        """
        Take a collapsed map and generate a file for histogramming in TOPCAT
        """

        assert(self.collapsed==True),'Run collapse first!'

        file='stack_pixels_%s_S.txt' % self.average
        f=open(file,'w')
        f.write('# i j stack_avg_muJy\n')
        for i in xrange(len(self.stack_avg[:])):
            for j in xrange(len(self.stack_avg[:,:])):
                f.write('%i %i %8.4f\n' % (i,j,Jy2muJy*self.stack_avg[i,j]))
        f.close()

        return


    def shuffle(self):
        """
        Reorder the stamps (to avoid SEx spatial correlations)
        """

        from numpy.random import shuffle
        shuffle(self.stamps)
        self.shuffle=True


    def runner(self):
        """
        Dump running mean and running std to pngs
        """

        import numpy
        import matplotlib.pyplot as plt
        from math import sqrt
        import os

        print 'NB squash.runner() produces a running MEAN only'
        print '   Make sure you ran squash.collapse() first...'

        assert(self.collapsed==True),'Run collapse first!'

        xcentre=self.x_max_S
        ycentre=self.y_max_S
        #online_avg=list(running_avg(self.stamps[:,xcentre,ycentre]))

        if self.simulate:
            file='running_mean_sim.txt'
        else:
            file='running_mean_dat.txt'
        f=open(file,'w')

        n,mean,M2=0,0,0
        running_means=[]; running_stds=[]
        for n,d in enumerate(self.stamps[:,xcentre,ycentre]):
            variance,n,mean,M2 = online_variance(d,n,mean,M2)
            online_std=Jy2muJy*numpy.sqrt(variance)/sqrt(n-1)
            online_mean=Jy2muJy*mean
            run='%i %f %f' % (n,online_mean,online_std)
            running_means.append(online_mean)
            running_stds.append(online_std)
            #print run
            f.write('%s\n' % run)
            if n % 10000 == 0: print '==> %6i <==' % n

        f.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(running_means)
        if self.simulate:
            pass
        else:
            ax.set_ylim((0.0,10.0))
        plt.title('mean_S/muJy (%s)' % os.path.basename(self.gals[0].filename))
        plt.xlabel('sample')
        plt.ylabel('mean_S/muJy')
        if self.simulate:
            plt.savefig('run_mean_sim_S.png')
        else:
            plt.savefig('run_mean_dat_S.png')
        plt.close()

        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.plot(running_stds)
        x=numpy.arange(1,10**5)
        y=numpy.arange(1,10**5)
        y=[20/sqrt(el) for el in y]
        plt.plot(x,y)
        plt.xscale('log')
        plt.yscale('log')
        plt.title('std_S/muJy (%s)' % os.path.basename(self.gals[0].filename))
        plt.xlabel('sample')
        plt.ylabel('std_S/muJy')
        if self.simulate:
            plt.savefig('run_std_sim_S.png')
        else:
            plt.savefig('run_std_dat_S.png')
        plt.close()

        return


    def plot(self,average=None):
        """
        Plot the median and mean maps to png
        """

        import matplotlib.pyplot as plt
        import os

        assert(self.collapsed==True),'Run collapse first!'

        #try:
        #    pass
        #except AttributeError:
        #    print 'First run cube.collapse()'
        #    return

        assert(self.average==average), 'Mismatch in averages! %s %s' \
            % (self.average,average)

        print '**Be careful to apply the same average to plot as for the collapse method!'

        if average==None:
            print 'Choice of average required (mean or median)!'
            return

        elif average=='mean':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cax=ax.imshow(self.stack_avg*Jy2muJy,interpolation='nearest')
            fig.colorbar(cax)
            plt.gca().invert_yaxis()
            plt.title('mean_S/muJy (%s)' % os.path.basename(self.gals[0].filename))
            plt.xlabel('RA/pix')
            plt.ylabel('Dec/pix')
            if self.simulate:
                plt.savefig('mean_sim_S.png')
            else:
                plt.savefig('mean_dat_S.png')
            plt.close()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            cax=ax.imshow(self.stack_avg_L,interpolation='nearest')
            fig.colorbar(cax)
            plt.gca().invert_yaxis()
            plt.title('mean_L/(W/Hz) (%s)' % os.path.basename(self.gals[0].filename))
            plt.xlabel('RA/pix')
            plt.ylabel('Dec/pix')
            if self.simulate:
                plt.savefig('mean_sim_L.png')
            else:
                plt.savefig('mean_dat_L.png')
            plt.close()

        elif average=='median':
            fig = plt.figure()
            ax = fig.add_subplot(111)
            cax=ax.imshow(self.stack_avg*Jy2muJy,interpolation='nearest')
            fig.colorbar(cax)
            plt.gca().invert_yaxis()
            plt.title('median_S/muJy (%s)' % os.path.basename(self.gals[0].filename))
            plt.xlabel('RA/pix')
            plt.ylabel('Dec/pix')
            if self.simulate:
                plt.savefig('median_sim_S.png')
            else:
                plt.savefig('median_dat_S.png')
            plt.close()

            fig = plt.figure()
            ax = fig.add_subplot(111)
            cax=ax.imshow(self.stack_avg_L,interpolation='nearest')
            fig.colorbar(cax)
            plt.gca().invert_yaxis()
            plt.title('median_L/(W/Hz) (%s)' % os.path.basename(self.gals[0].filename))
            plt.xlabel('RA/pix')
            plt.ylabel('Dec/pix')
            if self.simulate:
                plt.savefig('median_sim_L.png')
            else:
                plt.savefig('median_dat_L.png')
            plt.close()

        return

    def writeFITS(self,average=None):

        """
        Write map to FITS
        """

        import pyfits

        assert(self.collapsed==True),'Run collapse first!'

        assert(self.average==average), 'Mismatch in averages! %s %s' \
            % (self.average,average)

        print '**Be careful to apply the same average to plot as for the collapse method!'

        if average==None:
            print 'Choice of average required (mean or median)!'
            return

        if self.simulate:
            f1='sim'
        else:
            f1='dat'

        f='%s_%s_S.fits' % (average,f1)
        pyfits.writeto(f,self.stack_avg*Jy2muJy)

        f='%s_%s_L.fits' % (average,f1)
        pyfits.writeto(f,self.stack_avg_L*Jy2muJy)

        return

#-----------------------------------------------------------

class stamp(object):

    """
    """

    def __init__(self,gal,radio_map,radio_header,x_range,y_range,simulate,prng):

        import numpy,pyfits
        from math import cos,pi
        from numpy.random import RandomState

        self.gal=gal
        if simulate:
            #random.seed(1234) <-- this must only be called once!!
            x=prng.randint(x_range,radio_header['NAXIS1']-x_range)
            #print random.get_state()
            y=prng.randint(y_range,radio_header['NAXIS1']-y_range)
        else:
            x,y=WCSToPix(radio_header,gal.alpha,gal.delta)
        #print x,y,gal.alpha,gal.delta,cos(gal.delta*pi/180.0)

        xmin=iround(x-0.5*x_range)
        xmax=iround(x+0.5*x_range)
        ymin=iround(y-0.5*y_range)
        ymax=iround(y+0.5*y_range)

        #print xmin,x,xmax,ymin,y,ymax

        self.x=x
        self.y=y

        # -1's needed for 0 v 1 pixel-counting
        self.cutout = Jy2muJy*radio_map[0,0,ymin-1:ymax-1,xmin-1:xmax-1]
        self.mean   = Jy2muJy*numpy.average(self.cutout)
        self.median = Jy2muJy*numpy.median(self.cutout)
        self.std    = Jy2muJy*numpy.std(self.cutout)
        self.max    = Jy2muJy*numpy.max(self.cutout)
        self.sovern = self.max/self.std

        self.cutout_L = gal.factor*self.cutout
        self.mean_L   = Jy2muJy*numpy.average(self.cutout_L)
        self.median_L = Jy2muJy*numpy.median(self.cutout_L)
        self.std_L    = Jy2muJy*numpy.std(self.cutout_L)
        self.max_L    = Jy2muJy*numpy.max(self.cutout_L)
        self.sovern_L = self.max/self.std

        del radio_map,radio_header

        del gal,self.gal
        #self.plot() # <------- This was an expensive mistake!

        #hdu=pyfits.PrimaryHDU(self.cutout)
        #hdu.writeto('cutout2.fits')

        return

    def plot(self):

        import matplotlib.pyplot as plt

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.imshow(self.cutout[:,:],interpolation='nearest')
        plt.savefig('test.png')
        return

#-----------------------------------------------------------

def iround(x):
    """iround(number) -> integer
    Round a number to the nearest integer.
    from http://www.daniweb.com/software-development/python/threads/299459
    """
    return int(round(x) - .5) + (x > 0)

#-----------------------------------------------------------

def readMill(fmill='/Users/jtlz2/video/millennium/cats/zwart/mill-merged-jz-floc-not_star_basic_cutdown.fits'):
    """
    Read millennium simulation into memory
    cutdown columns are: ID, z, sfr, M_stellar
    """
    
    import pyfits

    fits = pyfits.open(fmill)
    simdata=fits[1].data
    simhead=fits[1].header
    fits.close()

    print 'Read %s' % fmill

    return simdata,simhead

#-----------------------------------------------------------

def readRadio(fradio='/Users/jtlz2/video/vla/VLA-VVDS-1.4GHZ.FITS'):
    """
    Read radio image into memory (whether map or noise map)
    """

    import pyfits

    #fradio='/Users/jtlz2/video/vla/VLA-VVDS-1.4GHZ.FITS'
    #fnoise='/Users/jtlz2/video/vla/backrms.fits'

    fits = pyfits.open(fradio)
    imdata=fits[0].data
    imhead=fits[0].header
    fits.close()

    print 'Read map %s' % fradio
    
    return imdata,imhead

#-----------------------------------------------------------

def writeRadio(radio_data,f='radio/radio_data.txt'):

    """
    Write radio image to ASCII pixel list
    """

    #f='radio_data.txt'
    d=open(f,'w')
    d.write('# m n pixval_muJy\n')

    for n,row in enumerate(radio_data[0,0,:,:]):
        row=radio_data[0,0,n,:]
        for m,pixel in enumerate(row):
            line='%4i %4i %f\n' % (m,n,Jy2muJy*pixel)
            d.write(line)
        d.flush()

    d.close()

    print 'Wrote %s' % f

    return

#-----------------------------------------------------------

def writeSNMap(radio_map,noise_map,fsovern='radio/radio_sovern.fits'):

    import pyfits
    import numpy

    sovern_map=radio_map/noise_map

    pyfits.writeto(fsovern,sovern_map)
    print 'Wrote %s' % fsovern

    return

#----------------------------------------------------------- 

def pixToWCS(imhead,x,y):
    """
    """
    alpha_centre = imhead['CRVAL1']
    delta_centre = imhead['CRVAL2']
    alpha_delta = imhead['CDELT1']
    delta_delta = imhead['CDELT2']
    x_centre = imhead['CRPIX1']
    y_centre = imhead['CRPIX2']
    n_steps_x = (x - x_centre)
    n_steps_y = (y - y_centre)
    alpha = alpha_centre + n_steps_x * alpha_delta
    delta = delta_centre + n_steps_y * delta_delta
    return alpha,delta

#----------------------------------------------------------- 

def WCSToPix(imhead,alpha,delta):
    """
    """
    from math import cos,pi

    print "***DO NOT USE THIS ROUTINE - USE astLib INSTEAD***"
    return None

    alpha_centre = imhead['CRVAL1']
    delta_centre = imhead['CRVAL2']
    alpha_delta = imhead['CDELT1']
    delta_delta = imhead['CDELT2']
    alpha_delta = alpha_delta / (cos(delta_centre*pi/180.0)) # cos(dec) term
    x_centre = iround(imhead['CRPIX1'])
    y_centre = iround(imhead['CRPIX2'])
    n_steps_x = iround((alpha - alpha_centre) / alpha_delta)
    n_steps_y = iround((delta - delta_centre) / delta_delta)
    x = x_centre + n_steps_x
    y = y_centre + n_steps_y

    #print x,y,alpha,delta

    return x,y

#----------------------------------------------------------- 

def plotRadio(imdata):

    import matplotlib.pyplot as plt
    from matplotlib import cm

    #dmn = 0.00001
    #dmx = 1000000
    #logspace = 10.**numpy.linspace(dmn, dmx, 100)
    #clevs = logspace

    fig = plt.figure()
    ax = fig.add_subplot(111)
    #plt.clim(1,5)

#    imvals = numpy.sort(imdata[0,0,:,:].flatten())
#    lo = imvals[0]
#    hi = imvals[-1]
#    steps = (imvals[::len(imvals)/256] - lo) / (hi - lo)
#    num_steps = float(len(steps))
#    interps = [(s, idx/num_steps, idx/num_steps) for idx, s in enumerate(steps)]
#    interps.append((1, 1, 1))
#    cdict = {'red' : interps,
#             'green' : interps,
#             'blue' : interps}
#    histeq_cmap = matplotlib.colors.LinearSegmentedColormap('HistEq', cdict)

    cax=ax.imshow(imdata[0,0,:,:])#,cmap=histeq_cmap)#cm.afmhot)#,clim=(0.3,1))#,vmin=0.005,vmax=0.020)#,norm=LogNorm(vmin=clevs[0], vmax=clevs[-1]))
    #cax.set_clim(0.005,0.020)
    #plt.afmhot()
    #cmap = mpl.cm.afmhot
    #norm = mpl.colors.Normalize(vmin=5, vmax=10)
    fig.colorbar(cax)
    #cax.set_clim(0,1)
    plt.title('VLA VVDS 1.4 GHz')
    plt.xlabel('RA/pix')
    plt.ylabel('Dec/pix')
    plt.savefig('radio.png')

    plt.close()


#-----------------------------------------------------------

def calculateVolume(z1,z2):

    """
    DEFUNCT, UNCHECKED & INCOMPLETE
    Calculate survey volume for Madau plot
    in Mpc**3
    """

    from cosmocalc import cosmocalc
    from math import pi

    V1=cosmocalc(z1)['VCM']
    V2=cosmocalc(z2)['VCM']
    dV=V2-V1

    degToSr=(180.0/pi)**2

    dV /= degToSr

    Avla = (2501 * 1.5/3600.0)**2

    V = Avla*dV
    
    return V

#-----------------------------------------------------------

def shufflePosns(gals):
    """
    Shuffle alpha and delta of galaxies
    Do this independently
    """

    print '****Shuffling RAs and Decs of sources independently'
    alphas=[gal.alpha for gal in gals]
    deltas=[gal.delta for gal in gals]
    numpy.random.shuffle(alphas)
    numpy.random.shuffle(deltas)
    for n,gal in enumerate(gals):
        gal.alpha=alphas[n]
        gal.delta=deltas[n]

    return gals

#----------------------------------------------------------- 

def scrambleMasses(gals,dm=None,seed=None):
    """
    Scramble the object masses by a gaussian of sigma dm
    """

    import numpy

    if seed is not None:
        numpy.random.seed(seed=seed)

    ngals=len(gals)
    noises=numpy.zeros(ngals)
    noises=numpy.random.normal(0.0,dm,ngals)

    n=0
    for n,gal in enumerate(gals):
        gal.mass_bon+=noises[n]
        if n > 0 and n % 10000 == 0: print '==> %6i <==' % n

    return gals

#----------------------------------------------------------- 

def scrambleRedshifts(gals,dz=None,seed=None):
    """
    Scramble the object redshifts by a gaussian of sigma dz
    """

    import numpy
    
    if seed is not None:
        numpy.random.seed(seed=seed)

    ngals=len(gals)
    noises=numpy.zeros(ngals)
    noises=numpy.random.normal(0.0,dz,ngals)

    factors=numpy.array([(1.0+gal.z_best) for gal in gals])
    noises*=factors

    n=0
    for n,gal in enumerate(gals):
        gal.z_best+=noises[n]
        if n > 0 and n % 10000 == 0: print '==> %6i <==' % n

    return gals

#----------------------------------------------------------- 

def updateCosmology(gals,lumdists):
    """
    Do the cosmological calculations
    One could make this a bound method of class galaxy (or something)
    """

    from cosmocalc import cosmocalc
    from math import pi,log10

    n=0

    #lumdists = [gal.lumdist for gal in gals]
    #z_bests = [gal.z_best for gal in gals]

 #   lumdists= map(cosmocalc(gal.z_best)['DL_Mpc'],\
  #                [gal.lumdist for gal in gals])

    #lumdists=calculateLumdist()

    for n,gal in enumerate(gals):
        #gal.lumdist = cosmocalc(gal.z_best)['DL_Mpc']
        gal.lumdist=queryLumdist(gal.z_best,lumdists)
        if gal.z_best < 0:
            gal.factor=0.0
        else:
            gal.factor = 9.52e18 * gal.lumdist * gal.lumdist \
                * 4 * pi * (1+gal.z_best)**(-0.2)
        #gal.factor = gal.lumdist * gal.lumdist \
        #    * (1+gal.z_best)**(-0.2)
        #gal.kabsmag=gal.kmag-5.0*((log10(gal.lumdist*1.0e6))-1.0)
        #gal.kabsmag=-1.0 # Dummy value for now
        #gal.age=cosmocalc(gal.z_best)['zage_Gyr']
        if n > 0 and n % 10000 == 0: print '==> %6i <==' % n

    #gal.factor=gal.factor * 9.52e18 * 4 * pi

    return gals

#-----------------------------------------------------------

def findNearest(array,value):
    """
    from http://stackoverflow.com/questions/2566412/find-nearest-value-in-numpy-array
    """
    import numpy
    idx=(numpy.abs(array-value)).argmin()
    #return array[idx]
    return idx

#-----------------------------------------------------------

def queryLumdist(z,cosmology):
    """
    """

    if z > cosmology[0,-1]:
        print 'redshift out of range - %f > %f' % (z,cosmology[0,-1])
        return None

    idx=findNearest(z,cosmology[0,:])
    lumdist=cosmology[1,idx]

    return lumdist

#-----------------------------------------------------------

def calculateLumdist(dump=False):
    """
    Checked against Ned Wright - works very well
    """

    from cosmocalc import cosmocalc
    import numpy

    zmax=10
    dz=0.001

    limit=int(float(zmax)/dz)
    array=numpy.zeros((2,limit))

    if dump:
        ff='cosmo/lumdist.txt'
        f=open(ff,'w')
    for n in xrange(limit):
        z=float(n)*dz
        array[0,n] = z
        array[1,n] = cosmocalc(z)['DL_Mpc']
        if dump: f.write('%f %f\n'%(array[0,n],array[1,n]))

    if dump:
        print '(Look in %s)' % ff
        f.close()

    print 'Generated %i z steps up to z=%3.2f' % (limit,array[0,-1])

    return array

#-----------------------------------------------------------

def weighted_avg_and_std(values, weights):
    """
    Returns the weighted average and standard deviation.
    values, weights -- Numpy ndarrays with the same shape.
    Based on http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    """
    import numpy
    from math import sqrt
    average = numpy.average(values, weights=weights)
    variance = numpy.dot(weights, (values-average)**2)/weights.sum()  # Fast and numerically precise
    return (average, sqrt(variance),weights.sum())

#----------------------------------------------------------- 

def realise(hals,imdata,imhead,imnoise,m=21000,n=100,root='realisations',
            noise_thresh=None,clip_thresh=None):
    """
    Generate n realisations of m randomly-positioned K-band sources
    Print to screen (grep for muJy)
    Write lots of pixel files..
    """

    import os
    import numpy
    from math import sqrt

    if not(os.path.exists(root)): os.mkdir(root)

    f='realisation_stats.txt'
    d=open(os.path.join(root,f),'w')
    head='# n seed mean_muJy std_muJy median_muJy\n'
    d.write(head)

    for realisn in range(n):
        seed=realisn
        e=sprint(hals[:m],seed=realisn)

        Ss,Ls,Ws,Ms,Ks,Kabs,stats=e.pixels(imdata,imhead,imnoise,simulate=True,\
                             noise_thresh=noise_thresh,clip_thresh=clip_thresh,\
                             directory=root)
        mean=numpy.average(Ss[2:])
        median=numpy.median(Ss[2:])
        std=numpy.std(Ss[2:])/sqrt(len(Ss[2:])-1)
        line='%i %i %f %f %f\n' % (realisn,seed,mean,std,median)
        d.write(line)
        if realisn > 0 and realisn % 10 == 0: d.flush()

    d.close()

    return

#-----------------------------------------------------------

def kbins(flag='all'):

    """
    Generate bins for S v K plots
    Only currently set up for 'all' [and bzk, eros classes]
    """

    if flag=='sbzk':
        bins=[19.0,20.0,21.0,21.6,22.0,22.6,23.0]
    elif flag=='pbzk':
        bins=[20.0,21.0,22.0,23.0]
    elif flag=='eros':
        bins=[19.0,20.0,20.5,21.0,21.5,22.0,22.5,23.0]
    else:
        bins=[17.5,18.5,19.0,19.5,20.0,20.5,21.0,21.5,22.0,22.5,23.0,23.5]

    return bins

#----------------------------------------------------------- 

def zSbins(flag='all'):

    """
    Generate bins for S v z plots
    """

    if flag=='sbzk':
        bins=[0.0,0.6,1.3,1.8,3.1,4.0] # 5.0 but we truncate at z=4.0
    elif flag=='pbzk':
        bins=[0.7,1.4,2.5,3.1,4.0] # 5.0 but we truncate at z=4.0
    else:
        return None

    return bins

#----------------------------------------------------------- 

def zLbins(flag='all'):

    """
    Generate bins for L v z plots
    """

    bins=[0.7,1.4,2.5] # 5.0 but we truncate at z=4.0

    return bins

#----------------------------------------------------------- 

def wrapKCut(gals,radio_map,radio_head,noise_map,kind='all',wrap_root='wrap',\
             noise_thresh=None,clip_thresh=None,masses='bon'):
    """
    Generate results for a range of (apparent) K bins
    """

    import os

    if not(os.path.exists(wrap_root)): os.mkdir(wrap_root)
    pixels_root=os.path.join(wrap_root,'pixels/')
    if not(os.path.exists(pixels_root)): os.mkdir(pixels_root)

    #bins=[x/2.0 for x in range(40,50)]
    bins=kbins(flag=kind)

    ff = os.path.join(wrap_root,'wrap_%s_K.dat' % kind)
    w=open(ff,'w')
    line='# nsrc K_mid S_1.4_muJy err_low_S_1.4_muJy err_high_S_1.4_muJy K_delta K_low K_high S_K_Jy\n'
    w.write(line)

    #noise_thresh=40.0 # muJy
    #clip_thresh=1.0e10 # sigma

    print '***Clipping noise at %4.2f muJy' % noise_thresh
    print '***Warning: do not clip sigma for the median (clip_thresh=%4.2f)' % clip_thresh

    for n in range(len(bins)-1):
        b1=bins[n]; b2=bins[n+1]
        print '*******************'
        print 'K bins %2.1f -> %2.1f' % (b2,b1)
        ggals=selectGals(gals,flag=kind,k_cut_faint=b2,k_cut_bright=b1,masses=masses)
        g=sprint(ggals)
        f='pixels_%s_K_%2.1f_%2.1f_noiseclip.dat' % (kind,b1,b2)
        f=os.path.join(pixels_root,f)
        Ss,Ls,Ws,Ms,Ks,Kabs,stats=g.pixels(radio_map,radio_head,noise_map,f=f,\
                       simulate=False,report=os.path.join(wrap_root,'wrap_%s_K_report.txt'%kind),\
                       noise_thresh=noise_thresh,\
                       clip_thresh=clip_thresh) # For median, don't clip sigma

        if stats['nsrc']==0:
            # When there are no sources in the given bin for some reason..
            print '***WARNING No sources in bin'
            print '  [Have you forgotten to set noise_thresh/clip_thresh??]'
            warn='***Skipping bin K %4.2f -> %4.2f' % (b1,b2)
            print warn
            continue

        K=(b1+b2)/2.0
        dK=(b2-b1)/2.0
        line='%6i %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %e \n' \
            % (stats['nsrc'],K,stats['S_median'],stats['S_median_err_low'],stats['S_median_err_high'],\
               dK,b1,b2,10**(-K/2.5))
        w.write(line)
        w.flush()

    w.close()

    return

#-----------------------------------------------------------

def wrapzCut(gals,radio_map,radio_head,noise_map,kind='all',wrap_root='wrap',\
             noise_thresh=None,clip_thresh=None,pixels_root='pixels'):
    """
    Generate results for a range of z bins
    (i.e., M is marginalized over)
    I THINK THIS IS NOW NOT USEFUL FOR ANY OF THE PLOTS
    IT WOULD'VE COVERED DUNNE FIG. 12, BUT THAT CAN NOW
        BE HANDLED BY wrapMzCut
    """

    import os

    if kind == 'pbzk' or kind == 'sbzk':
        pass
    else:
        print '%s not supported!' % kind
        return

    bins=zSbins(flag=kind)

    ff = os.path.join(wrap_root,'wrap_%s_z.dat' % kind)
    w=open(ff,'w')
    line='# z_mid S_1.4_muJy err_S_1.4_muJy L_WperHz err_L_WperHz z_delta z_low z_high M_st_peryr err_M_st_peryr\n'
    w.write(line)

    for n in range(len(bins)-1):
        b1=bins[n]; b2=bins[n+1]
        print '*******************'
        print 'z bins %2.1f -> %2.1f' % (b1,b2)
        ggals=selectGals(gals,flag=kind,z_cut_low=b1,z_cut_high=b2,k_cut_faint=23.5)
        g=sprint(ggals)
        f='pixels_%s_z_%2.1f_%2.1f_noiseclip.dat' % (kind,b1,b2)
        f=os.path.join(pixels_root,f)
        Ss,Ls,Ws,Ms,Ks,Kabs,stats=g.pixels(radio_map,radio_head,noise_map,f=f,\
                       simulate=False,\
                       report=os.path.join('wrap','wrap_%s_z_report.txt'%kind),\
                       noise_thresh=noise_thresh,\
                       clip_thresh=clip_thresh) # For median, don't clip sigma


        line='%2.1f %2.1f %2.1f %e %e %5.3f %5.3f %5.3f %5.3f %5.3f\n' \
            % ((b1+b2)/2.0,stats['S_median'],stats['S_median_err'],\
               stats['L_median'],stats['L_median_err'],\
               (b2-b1)/2.0,b1,b2,\
               stats['M_median'],stats['M_median_err'])
        w.write(line)
        w.flush()

    w.close()

    return

#-----------------------------------------------------------

def zipMill(mdata):

    """
    Carry out the slow zip/transpose steps on the Millennium data
    """

    import numpy

    data=numpy.array(mdata)
    zipped=map(list,zip(*data))
    e=numpy.array(zipped)
    transposed=numpy.transpose(e)
    del data,zipped,e

    return transposed


#-----------------------------------------------------------

def generateMillGrid(zbins=[],Mbins=[],mdata=None,zippy=None):

    """
    Take millennium simulation data and grid it into provided M,z bins
    stacker.generateMillGrid(zbins=[0.2,0.7,1.2,1.6,2.3,2.8,3.9],
        Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0])
    """

    import numpy
    from math import log10

    mfile='/Users/jtlz2/video/millennium/cats/zwart/mill-merged-jz-floc-not_star_basic_cutdown.fits'
    mfile='/Users/jtlz2/video/millennium/cats/zwart/K23.5/mill-merged-jz-clean-plus-cutdown.fits'
    if mdata is None:
        mdata,mhead=readMill(mfile)

    if zippy is None:
        full=zipMill(mdata)
    else:
        full=zippy

    #(idcol,zcol,sfrcol,Mcol)=range(4)
    #(zcol,kcol,Mcol,sfrcol)=range(4)
    (zcol,ssfrcol,Mcol,sfrcol)=range(4)
    print '# nsrc zc zlow zhigh Mc Mlow Mhigh zmed Mmed sfr err_sfr ssfr err_ssfr log_ssfr log_ssfr_err'
    for n in range(len(zbins)-1):
        zb1=zbins[n]; zb2=zbins[n+1]
        for j in range(len(Mbins)-1):
            Mb1=Mbins[j]; Mb2=Mbins[j+1]
            f=full[full[:,zcol]>zb1]
            #print zb1,zb2,Mb1,Mb2,zcol,Mcol,numpy.shape(f),numpy.shape(full)
            try:
                assert(numpy.size(f)>0)
            except AssertionError:
                print '# ***WARNING(1): Empty bin %5.3f %5.3f' % ((zb1+zb2)/2.0,(Mb1+Mb2)/2.0)
                continue

            g=f[f[:,zcol]<zb2]
            i=g[g[:,Mcol]>Mb1]
            #gavo.mpa-garching.mpg.de/Millennium/Help?page=databases/henriques2012a/database:
            # Ignore masses < 1e9 (i.e. log10 M < 0.1)
            #k=i[i[:,Mcol]>0.10] # Condition from Henriques 2012a note on Mill. website
            try:
                assert(numpy.size(i)>0)
            except AssertionError:
                print '# ***WARNING(2): Empty bin %5.3f %5.3f' % ((zb1+zb2)/2.0,(Mb1+Mb2)/2.0)
                continue
            k=i[i[:,Mcol]>9.0] # Condition from Henriques 2012a note on Mill. website
            #except IndexError:

            binned=k[k[:,Mcol]<Mb2]
            del f,g,i,k

            h=1.0 # <--- Andreas says h=1.0 corresponds to H_0=73 (Mill. sims) [so ignore]
            # units of M_stellar are 10^10 Msun /h
            ssfr_arr=1.0e9*binned[:,sfrcol]/(numpy.power(10,binned[:,Mcol])) # /Gyr
            ssfr2_arr=binned[:,ssfrcol]

            ssfr=numpy.median(ssfr_arr)   # }
            ssfr2=numpy.median(ssfr2_arr) # } Confirmed these two come out the same

            err_ssfr=err(ssfr_arr)[2]

            sfr_arr=binned[:,sfrcol]
            sfr=numpy.median(sfr_arr)
            err_sfr=err(sfr_arr)[2]

            if ssfr > 0.0:
                log_ssfr = log10(ssfr)
            else:
                log_ssfr = float('nan')
            log_ssfr_err = log10(euler)*err_ssfr/ssfr

            #assert(numpy.median(ssfr)==numpy.median(ssfr2)), \
            #    'ssfr discrepancy %f %f' % (numpy.median(ssfr),numpy.median(ssfr2))
            nsrc=numpy.size(ssfr_arr)
            zmed=numpy.median(binned[:,zcol])
            Mmed=numpy.median(binned[:,Mcol])
            zc=(zb1+zb2)/2.0
            Mc=(Mb1+Mb2)/2.0
            if nsrc <= 2:
                line = '# *** WARNING: nsrc < %i (%5.3f,%5.3f)' \
                    % (nsrc,zc,Mc)
            else:
                line ='%8i %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %e %e %e %e' \
                    % (nsrc,zc,zb1,zb2,Mc,Mb1,Mb2,zmed,Mmed,sfr,\
                       err_sfr,ssfr,err_ssfr,log_ssfr,log_ssfr_err)
            print line
            del binned,ssfr_arr,sfr_arr

    return


#-----------------------------------------------------------

def wrapMzCut(gals,radio_map,radio_head,noise_map,kind='all',wrap_root='test',\
             noise_thresh=None,clip_thresh=None,pixels_root='test/pixels_Mz',\
             include_detections=False,masses='bon',zbins=[],Mbins=[],simulate=False,\
             KLIM=23.5,kabsbins=[-100.0,0.0],upsampling_factor=None,\
             upsampled_map=None,upsampled_noise=None,use_radio_posn=False,\
             do_sim_posns=None,scramble_posns=None,noise_only=None,\
             scramble_posns_seed=None,scramble_posns_mask=None):
    """
    Generate results for a range of M and z bins
    This routine can be used for Dunne Figs 12 and 13 - just vary the bins (below)
    [Plot these using gnuplot scripts]
    ALSO now handles Kabs bins
    pixels_root is now ignored if supplied - is generated from wrap_root
    """

    import os


    report=os.path.join(wrap_root,'wrap_%s_z_M_report.txt'%kind)

    pixels_root=os.path.join(wrap_root,'pixels_Mz')

    if not(os.path.exists(wrap_root)): os.mkdir(wrap_root)
    if not(os.path.exists(pixels_root)): os.mkdir(pixels_root)
    if os.path.exists(os.path.join(report)): os.remove(report)

    #
    #*** ->   zbins=[0.0,0.5,1.0,1.5,2.0,2.5,3.5]#,3.5,4.0]
    #zbins=[0.25,0.75,1.25,1.75,2.25,2.75,3.25]#,3.5,4.0]
    # zbins=[0.2,0.725,1.175,1.625,2.175,2.9] #These are only approximations to Dunne's bins
    #zbins=[0.0,1.0,10.0]
    #zbins=[0.0,4.0]   # ***THE UPPER LIMIT OF 4.0 IS IMPORTANT
    #zbins=[-1.0,1.4,4.0]
    #*** ->   zbins=[0.2,0.7,1.2,1.6,2.3,2.8,3.9]

    # Strictly - should have overlapping bins to match Dunne - skip for now
    #*** ->   Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    #Mbins=[9.0,12.0]
    #Mbins=[9.0,10.5,12.0]
    #Mbins=[-90.0,100.0]
    #Mbins=[7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    #Mbins=[5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]

    #kabsbins=[-100.0,0.0]
    #kabsbins=[-26.0,-25.0,-24.0,-23.0,-22.0,-21.0,-20.0]
    #kabsbins=[-100,-25.0]

    ff = os.path.join(wrap_root,'wrap_%s_Mz.dat' % kind)
    w=open(ff,'w')
    line= '%s%s%s%s%s%s%s%s%s'%('# nsrc z_mid S_1.4_muJy err_S_1.4_muJy L_WperHz err_L_WperHz ', \
                 'z_delta z_low z_high M_st_peryr err_M_st_peryr M_mid M_delta ', \
                 'M_low M_high SFR_c_msolperyr err_SFR_c_msolperyr K K_err ', \
                 'Kabs Kabs_err SFR_y_msolperyr err_SFR_y_msolperyr SSFR_c_peryr ', \
                 'err_SSFR_c_peryr SSFR_y_peryr err_SSFR_y_peryr z err_z ', \
                 'log_SSFR_c_peryr err_log_SSFR_c_peryr ', \
                 'log_SSFR_y_peryr err_log_SSFR_y_peryr ', \
                 'SFR_b SSFR_b SFR_g SSFR_g err_SSFR_b err_SSFR_g log_SSFR_b log_SSFR_g ', \
                 'err_log_SSFR_b err_log_SSFR_g err_SFR_b err_SFR_g\n')
    w.write(line)

    for n in range(len(zbins)-1):
        zb1=zbins[n]; zb2=zbins[n+1]

        for j in range(len(Mbins)-1):
            Mb1=Mbins[j]; Mb2=Mbins[j+1]

            for m in range(len(kabsbins)-1):
                Kb1=kabsbins[m]; Kb2=kabsbins[m+1]

                print '*******************'
                print 'z bins %2.1f -> %2.1f, M bins %2.1f -> %2.1f, Kabs bins %2.1f -> %2.1f'\
                    % (zb1,zb2,Mb1,Mb2,Kb1,Kb2)
                ggals=selectGals(gals,flag=kind,z_cut_low=zb1,\
                                 z_cut_high=zb2,k_cut_faint=KLIM,\
                                 M_stellar_low=Mb1,M_stellar_high=Mb2,\
                                 kabs_cut_faint=Kb2,kabs_cut_bright=Kb1,\
                                 include_detections=include_detections,\
                                 masses=masses)

                                 #for gal in ggals:
                                 #if gal.z_best >= 4.0: print gal.z_best

                print 'n in bin %i' % int(len(ggals))
                #ggals=ggals[:int(len(ggals)/2.0)]
                #print len(ggals)
                g=sprint(ggals)
                f='pixels_%s_z_%2.1f_%2.1f_Ms_%2.1f_%2.1f_Kabs_%2.1f_%2.1f_noiseclip.dat' \
                    % (kind,zb1,zb2,Mb1,Mb2,Kb1,Kb2)
                f=os.path.join(pixels_root,f)

                Ss,Ls,Ws,Ms,Ks,Kabs,stats=g.pixels(radio_map,radio_head,noise_map,f=f,\
                                       simulate=simulate,\
                                       report=report,\
                                       noise_thresh=noise_thresh,\
                                       clip_thresh=clip_thresh,\
                                       masses=masses,
                                       upsampling_factor=upsampling_factor,\
                                       upsampled_map=upsampled_map,
                                       upsampled_noise=upsampled_noise,\
                                       use_radio_posn=use_radio_posn,\
                                       do_sim_posns=do_sim_posns,\
                                       scramble_posns=scramble_posns,\
                                       noise_only=noise_only,\
                                       scramble_posns_seed=scramble_posns_seed,\
                                       scramble_posns_mask=scramble_posns_mask) # For median, don't clip sigma
                if stats['nsrc']==0:
                    # When there are no sources in the given bin for some reason..
                    print '***WARNING No sources in bin'
                    print '  [Have you forgotten to set noise_thresh/clip_thresh??]'
                    warn='***Skipping bin z %4.2f -> %4.2f : M %4.2f -> %4.2f : Kabs %4.2f -> %4.2f' \
                        % (zb1,zb2,Mb1,Mb2,Kb1,Kb2)
                    print warn
                    #line='# %s z %4.2f %4.2f  M %4.2f %4.2f %4.2f %4.2f\n' \
                    #    % ('0?',zb1,zb2,Mb1,Mb2,Kb1,Kb2)
                    #w.write(line)
                    continue

                if not simulate:
                    line='%6i %5.3f %5.3f %5.3f %e %e %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %e %e %e %e %5.3f %5.3f %e %e %e %e %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %e %e %e %e %e %e %e %e %e %e %e %e\n' \
                        % (stats['nsrc'],(zb1+zb2)/2.0,
                       stats['S_median'],stats['S_median_err'],
                       stats['L_median'],stats['L_median_err'],
                       (zb2-zb1)/2.0,zb1,zb2,
                       stats['M_median'],stats['M_median_err'],
                       (Mb1+Mb2)/2.0,(Mb2-Mb1)/2.0,Mb1,Mb2,
                       stats['SFR_condon_median'],stats['SFR_condon_median_err'],
                       stats['K_median'],stats['K_median_err'],
                       stats['Kabs_median'],stats['Kabs_median_err'],
                       stats['SFR_yun_median'],stats['SFR_yun_median_err'],
                       stats['SSFR_condon_median'],stats['SSFR_condon_median_err'],
                       stats['SSFR_yun_median'],stats['SSFR_yun_median_err'],
                       stats['z_median'],stats['z_median_err'],
                       stats['log_SSFR_condon_median'],stats['log_SSFR_condon_median_err'],\
                       stats['log_SSFR_yun_median'],stats['log_SSFR_yun_median_err'], \
                       stats['SFR_bon_median'],stats['SSFR_bon_median'], \
                       stats['SFR_gio_median'],stats['SSFR_gio_median'], \
                       stats['SSFR_bon_median_err'],stats['SSFR_gio_median_err'], \
                       stats['log_SSFR_bon_median'],stats['log_SSFR_gio_median'], \
                       stats['log_SSFR_bon_median_err'],stats['log_SSFR_gio_median_err'], \
                       stats['SFR_bon_median_err'],stats['SFR_gio_median_err'])

                w.write(line)
                w.flush()

    w.close()

    if stats['nsrc'] > 0 and not simulate:
        #print '%s  &  %5i  & %4.2f  & $\pm$ &  %4.2f $\pm$ %4.2f  &  %4.2e $\pm$ %4.2e  &  %4.2f $\pm$ %4.2f  &  %4.2f $\pm$ %4.2f  &  %4.2f \\\\ \n' % \
        #    (kind,stats['nsrc'],stats['frac_detected'],stats['S_median'],\
        #     stats['S_median_err'],stats['L_median'],stats['L_median_err'],\
        #     stats['SFR_condon_median'],stats['SFR_condon_median_err'],\
        #     stats['SSFR_condon_median'],stats['SSFR_condon_median_err'],stats['z_median'])

        # Symmetric error bars
        if stats['S_median_err_low'] == stats['S_median_err_high']:
            assert(stats['S_median_err_low']==stats['S_median_err_high']\
                   ==stats['S_median_err']), 'Problem with error bars!'
            print '%s  &  %5i  & %5.3f  & $\pm$ &  $%5.3f\pm %5.3f$ &  $%5.3e\pm %5.3e$  &  $%5.3f\pm %5.3f$  &  $%5.3f\pm %5.3f$  &  %5.3f \\\\ \n' % \
                (kind,stats['nsrc'],100.0*stats['frac_detected'],\
                stats['S_median'],stats['S_median_err'],\
                stats['L_median'],stats['L_median_err'],\
                stats['SFR_condon_median'],stats['SFR_condon_median_err'],\
                stats['SSFR_condon_median'],stats['SSFR_condon_median_err'],\
                stats['z_median'])

        # Asymmetric error bars
        else:
            print '%s  &  %5i  & %5.3f  & $\pm$ &  $%5.3f_{-%5.3f}^{+%5.3f}$ &  $%5.3e_{-%5.3e}^{+%5.3e}$  &  $%5.3f_{-%5.3f}^{+%5.3f}$  &  $%5.3f_{-%5.3f}^{+%5.3f}$  &  %5.3f \\\\ \n' % \
                (kind,stats['nsrc'],stats['frac_detected'],stats['S_median'],\
                stats['S_median_err_low'],stats['S_median_err_high'],\
                stats['L_median'],\
                stats['L_median_err_low'],stats['L_median_err_high'],\
                stats['SFR_condon_median'],\
                stats['SFR_condon_median_err_low'],stats['SFR_condon_median_err_high'],\
                stats['SSFR_condon_median'],\
                stats['SSFR_condon_median_err_low'],stats['SSFR_condon_median_err_high'],\
                stats['z_median'])


    return

#-----------------------------------------------------------


if __name__ == "__main__":

    import sys

    import stacker
    reload(stacker)
    # Read the radio data
    imdata,imhead=stacker.readRadio(fradio='/Users/jtlz2/video/vla/VLA-VVDS-1.4GHZ.FITS')
    # Read the noise map
    noisedata,noisehead=stacker.readRadio(fradio='/Users/jtlz2/video/vla/backrms.fits')
    # Read the K-band catalogue (up to z=5)
    gals=stacker.readCatalogue(do_shuffle=False,version=3)

    #c=stacker.cube(gals,imdata,imhead,simulate)
    #d=stacker.squash(c)
    #d.collapse(average='mean')
    #d.shuffle()
    #d.runner()
    #d.plot(average='mean')
    #d.writeFITS(average='mean')

    # Set up the D_L curve (do once)
    cosmo=stacker.calculateLumdist()

    # **** December 2013 runs:
    # Now do these next two steps below:
    
    # Optionally scramble the redshifts
    gals=stacker.scrambleRedshifts(gals,dz=0.13,seed=1234)

    # Use the D_L curve to generate the Dunne cosmology factors
    gals=stacker.updateCosmology(gals,cosmo)

    # Skip to ** below

    # Select galaxy type here
    #    hals=stacker.selectGals(gals,flag=None)

    # Instantiate pixel-stack object(?!)
    #       and setup the simulation seed (even if not required)
    #    h=stacker.sprint(hals)
    #    Ss,Ls,Ws,Ms,Ks,Kabs,stats=h.pixels(imdata,imhead,noisedata,\
    #        f='pixels.dat',simulate=False)

    # Wrapper for generating data for plots
    # postprocess using gnuplot scripts
    stacker.wrapKCut(gals,imdata,imhead,noisedata,kind='eros')
    stacker.wrapzCut(gals,imdata,imhead,noisedata,kind='pbzk')

    stacker.wrapKCut(gals,imdata,imhead,noisedata,kind='not_star',wrap_root='production2_c',noise_thresh=40,clip_thresh=1.0e10)
    stacker.wrapKCut(gals,imdata,imhead,noisedata,kind='sbzk',wrap_root='production2_d',noise_thresh=40,clip_thresh=1.0e10)
    stacker.wrapKCut(gals,imdata,imhead,noisedata,kind='pbzk',wrap_root='production2_e',noise_thresh=40,clip_thresh=1.0e10)
    stacker.wrapKCut(gals,imdata,imhead,noisedata,kind='eros',wrap_root='production2_f',noise_thresh=40,clip_thresh=1.0e10)

    # Generate SvK plot
    mass='bon'
    stacker.wrapKCut(gals,imdata,imhead,noisedata,kind='not_star',wrap_root='production3_z',noise_thresh=40,clip_thresh=1.0e10,masses=mass)


    # For the median, set clip_thresh way high
    # For the mean, set clip_thresh to 5.0 sigma

    # ** Generate Mz grid (and everything else..) for plotting
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,kind='not_star',\
                      noise_thresh=40.0,clip_thresh=1.0e10,\
                      wrap_root='SSFRvMz',pixels_root='SSFRvMz/pixels_Mz',\
                      include_detections=False,masses='bon')

    # Compare with jarvis/4JZ_K23_alltogether.dat
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='SSFRvMz_test2d_all4b',\
        pixels_root='SSFRvMz_test2d_all4b/pixels_Mz',zbins=[0,4],\
        Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0])

    # Count the objects for each gzK-classification
    # --->> CHECK the Mbins range is OK
    mass='bon'
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='production1_o',\
        pixels_root='production1_o/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='not_star',wrap_root='production1_x',\
        pixels_root='production1_x/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='production2_g',\
        pixels_root='production2_g/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass,include_detections=True)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sgzk',wrap_root='production1_p',\
        pixels_root='production1_p/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='pgzk',wrap_root='production1_q',\
        pixels_root='production1_q/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='pgzk',wrap_root='production1_y',\
        pixels_root='production1_y/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ngzk',wrap_root='production1_r',\
        pixels_root='production1_r/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='ngzk',wrap_root='production1_z',\
        pixels_root='production1_z/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)


    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='sgzk',wrap_root='production1_s',\
        pixels_root='production1_s/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    # Check error bars by halving the number of sources
    # This cut actually fails (why??), so instead hack inside wrapMzCut

    #    half_gals=gals[:int(len(gals)/2.0)]
    #half_gals=gals
    #print len(gals), len(half_gals)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='production2_h',\
        pixels_root='production2_h/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='not_star',wrap_root='production2_i',\
        pixels_root='production2_i/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    # ================================================================================
    # Produce Dunne Table 2 again - counts and SSFR results
    mass='bon'
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='production2_q',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='not_star',wrap_root='production2_r',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sgzk',wrap_root='production2_s',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='sgzk',wrap_root='production2_t',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='pgzk',wrap_root='production2_u',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='pgzk',wrap_root='production2_v',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ngzk',wrap_root='production2_w',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='ngzk',wrap_root='production2_x',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='eros',wrap_root='production2_y',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='eros',wrap_root='production2_z',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    # The next two (3_a,3_b) have the wrong flags - corrected in time for for 3_t,3_u
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='drgs',wrap_root='production3_a',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='drgs',wrap_root='production3_b',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='star',wrap_root='production3_c',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='star',wrap_root='production3_d',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='drgs',wrap_root='production3_t',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='drgs',wrap_root='production3_u',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    # Asymmetric error bars
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='drgs',wrap_root='production3_v',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='drgs',wrap_root='production3_w',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='production3_x',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='not_star',wrap_root='production3_y',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    # ================================================================================
    # Splits by MOD_BEST

    mass='bon'
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production3_e',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='ell',wrap_root='production3_f',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbc',wrap_root='production3_m',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='sbc',wrap_root='production3_n',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='scd',wrap_root='production3_g',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='scd',wrap_root='production3_h',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='irr',wrap_root='production3_i',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='irr',wrap_root='production3_j',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production3_k',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='sbn',wrap_root='production3_l',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass)

    mass='bon'
    zbins=[0.2,0.7,1.2,1.6,2.3,2.8,3.9]
    Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production3_o',\
                      zbins=zbins,Mbins=Mbins,masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbc',wrap_root='production3_p',\
                      zbins=zbins,Mbins=Mbins,masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='scd',wrap_root='production3_q',\
                      zbins=zbins,Mbins=Mbins,masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='irr',wrap_root='production3_r',\
                      zbins=zbins,Mbins=Mbins,masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production3_s',\
                      zbins=zbins,Mbins=Mbins,masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='drgs',wrap_root='production3_y',\
                      zbins=zbins,Mbins=Mbins,masses=mass)


    # [120511] Final runs
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='production4_test',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production4_b',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='ell',wrap_root='production4_c',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='mid',wrap_root='production4_d',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='mid',wrap_root='production4_e',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production4_f',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='sbn',wrap_root='production4_g',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production4_h',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='all',wrap_root='production4_i',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='star',wrap_root='production4_j',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='star',wrap_root='production4_k',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')

    # TEST OF PIXEL OFFSETS
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='all',wrap_root='production4_t',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production4_u',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')

    # Search for pixels that fail the 40-muJy noise cut
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=1.0e10,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production4_v',\
                      zbins=[0,5],Mbins=[-90.0,100.0],masses='bon')


    # Generate SvK plot
    stacker.wrapKCut(gals,imdata,imhead,noisedata,kind='all',wrap_root='production4_l',\
                     noise_thresh=40,clip_thresh=1.0e10,masses='bon')

    # Generate SSFR v M,z plots
    #zbins=[0.2,0.7,1.2,1.6,2.3,2.8,3.9]
    zbins=[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
    Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production4_m',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='mid',wrap_root='production4_n',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production4_o',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production4_q',\
                      zbins=zbins,Mbins=Mbins,masses='bon')

    # Generate Mz grid for millennium simulation data -> 4_a [Fixed bins]
    # Note that resolution of simulation -> Mbin_0 = 9.0
    mfile='/Users/jtlz2/video/millennium/cats/zwart/K23.5/mill-merged-jz-clean-plus-cutdown.fits'
    mdata,mhead=stacker.readMill(mfile)
    zippy=stacker.zipMill(mdata)
    stacker.generateMillGrid(zbins=[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0],
        Mbins=[9.0,9.5,10.0,10.5,11.0,11.5,12.0],zippy=zippy)

    # Generate Mz grid for millennium simulation data -> 4_r [Dunne-style bins]
    stacker.generateMillGrid(zbins=[0.2,0.7,1.2,1.6,2.3,2.8,3.9],
        Mbins=[9.0,9.5,10.0,10.5,11.0,11.5,12.0],zippy=zippy)

    # Generate Mz single value for millennium simulation data -> 4_s
    stacker.generateMillGrid(zbins=[0,5],Mbins=[9.0,100.0],zippy=zippy)

    # [Dunne bins - I think these look better]
    zbins=[0.2,0.7,1.2,1.6,2.3,2.8,3.9]
    Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production4_p',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    # [/end Dunne bins/]


    # ================================================================================
    # Test out simulation of random positions
    mass='bon'
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='not_star',wrap_root='production2_test',\
                      zbins=[0,4],Mbins=[-90.0,100.0],masses=mass,simulate=True)


    mass='gio'
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sgzk',wrap_root='production1_t',\
        pixels_root='production1_t/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='sgzk',wrap_root='production1_u',\
        pixels_root='production1_u/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    mass='bon'
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='eros',wrap_root='production1_v',\
        pixels_root='production1_v/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='eros',wrap_root='production1_w',\
        pixels_root='production1_w/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='drgs',wrap_root='production2_a',\
        pixels_root='production2_a/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=5.0,kind='drgs',wrap_root='production2_b',\
        pixels_root='production2_b/pixels_Mz',zbins=[0,4],Mbins=[-90.0,100.0],\
        masses=mass)

    mass='bon'
    zbins=[0.2,0.7,1.2,1.6,2.3,2.8,3.9]
    Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='not_star',wrap_root='production2_j',\
        pixels_root='production2_j/pixels_Mz',zbins=zbins,Mbins=Mbins,\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sgzk',wrap_root='production2_k',\
        pixels_root='production2_k/pixels_Mz',zbins=zbins,Mbins=Mbins,\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='pgzk',wrap_root='production2_l',\
        pixels_root='production2_l/pixels_Mz',zbins=zbins,Mbins=Mbins,\
        masses=mass)
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ngzk',wrap_root='production2_m',\
        pixels_root='production2_m/pixels_Mz',zbins=zbins,Mbins=Mbins,\
        masses=mass)

    zbins=[0,4]
    Mbins=[-90.0,100.0]


    #SIMULATION:
    # Simulate n realisations of randomly-positioned K-band sources
    stacker.realise(hals,imdata,imhead,noisedata,n=100,m=21000,\
                    noise_thresh=40.0,clip_thresh=1.0e10)
    stacker.realise(gals,imdata,imhead,noisedata,n=100,m=1000,\
                    noise_thresh=40.0,clip_thresh=1.0e10,root='test_rx2')
    stacker.realise(gals,imdata,imhead,noisedata,n=40,m=1000,\
                    noise_thresh=40.0,clip_thresh=1.0e10,root='test_rx')
    stacker.realise(gals,imdata,imhead,noisedata,n=100,m=10000,\
                    noise_thresh=40.0,clip_thresh=1.0e10,root='test_rx3')
    stacker.realise(gals,imdata,imhead,noisedata,n=100,m=40000,\
                    noise_thresh=40.0,clip_thresh=1.0e10,root='test_rx4')

    # Pretty final sims (1000 realns takes about 8 hours - cut short before then):
    stacker.realise(gals,imdata,imhead,noisedata,n=1000,m=40000,\
                    noise_thresh=40.0,clip_thresh=1.0e10,root='test_rx5')
    stacker.realise(gals,imdata,imhead,noisedata,n=1000,m=40000,\
                    noise_thresh=40.0,clip_thresh=5.0,root='test_rx6')


                    #    log10(SFR_c_msolperyr/pow(10,M_mid))
                    #    = log10(SFR_c_msolperyr) - log10(10^M_mid)
                    #    = log10(SFR_c_msolperyr) - M_mid

    # ================================================================================
    # December 2013 runs after referee's comments

    imdata,imhead=stacker.readRadio(fradio='/Users/jtlz2/video/vla/VLA-VVDS-1.4GHZ.FITS')
    # Read the noise map
    noisedata,noisehead=stacker.readRadio(fradio='/Users/jtlz2/video/vla/backrms.fits')
    # Read the K-band catalogue
    gals=stacker.readCatalogue(do_shuffle=False,version=3)

    # Set up the D_L curve (do once)
    cosmo=stacker.calculateLumdist()

    # Use the D_L curve to generate the Dunne cosmology factors
    gals=stacker.updateCosmology(gals,cosmo)

    # Median runs only - scrambled redshifts

    # 0 < z < 3 --- 9.0 < M < 11.5 --- cf 4_h etc.
    # all 49604 <z>=1.348
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production5_h',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')
    # ell 9900 <z>=1.442
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production5_i',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')
    # mid 33747 <z>=1.334
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='mid',wrap_root='production5_j',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')
    # sbn 5957 <z>=1.207
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production5_k',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')


    # cf 5_m 5_n 5_o 5_q
    #zbins=[0.0,0.5,1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0]
    #Mbins=[6.5,7.0,7.5,8.0,8.5,9.0,9.5,10.0,10.5,11.0,11.5,12.0]
    zbins=[0.0,0.5,1.0,1.5,2.0,2.5,3.0]
    Mbins=[9.0,9.5,10.0,10.5,11.0,11.5]
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production5_q',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production5_m',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='mid',wrap_root='production5_n',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(gals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production5_o',\
                      zbins=zbins,Mbins=Mbins,masses='bon')

    # Scramble the redshifts
    zgals=stacker.scrambleRedshifts(gals,dz=0.13,seed=1234)

    # Use the D_L curve to generate the Dunne cosmology factors
    zgals=stacker.updateCosmology(zgals,cosmo)

    # 0 < z < 3 --- 9.0 < M < 11.5 --- cf 4_h etc.
    # all 48844 <z>=1.310
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production5_r',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')
    # ell 9766 <z>=1.397
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production5_s',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')
    # mid 33189 <z>=1.297
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='mid',wrap_root='production5_t',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')
    # sbn 5889 <z>=1.229
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production5_u',\
                      zbins=[0,3],Mbins=[9.0,11.5],masses='bon')

    zbins=[0.0,0.5,1.0,1.5,2.0,2.5,3.0]
    Mbins=[9.0,9.5,10.0,10.5,11.0,11.5]
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production5_v',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='ell',wrap_root='production5_w',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='mid',wrap_root='production5_x',\
                      zbins=zbins,Mbins=Mbins,masses='bon')
    stacker.wrapMzCut(zgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='sbn',wrap_root='production5_y',\
                      zbins=zbins,Mbins=Mbins,masses='bon')

    # Scramble the redshifts AND masses
    zgals=stacker.scrambleRedshifts(gals,dz=0.13,seed=1234)
    zgals=stacker.updateCosmology(zgals,cosmo)
    mgals=stacker.scrambleMasses(zgals,dm=0.1,seed=1235)
    stacker.wrapMzCut(mgals,imdata,imhead,noisedata,noise_thresh=40,\
                      clip_thresh=1.0e10,kind='all',wrap_root='production5_z',\
                      zbins=zbins,Mbins=Mbins,masses='bon')

    sys.exit(0)

