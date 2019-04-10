#!/usr/bin/env python

"""
Run as

./reconstruct.py chains_140107a

"""

import os,sys
import importlib
import numpy
from scipy import stats
from utils import sqDeg2sr,calculate_confidence,peak_confidence,medianArray,gaussian
import pymultinest
import lumfunc # This leads to lumfunc.data (not used here)
from math import log10,pi,e
import matplotlib.pyplot as plt
import pylab
import pofd

param_file=sys.argv[-1]
setf='%s.settings' % param_file
#try:
#    execfile(param_file)
#except IOError:
#    from settings import *
#setf='chains_140108a.settings'
#from chains_140108a.settings import *

#-------------------------------------------------------------------------------

def main():

    """
    """

    global datafile,NSKADS_RESCALING

    # Import the settings variables
    print 'Settings file is %s' % param_file

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    # Some global plot settings
    lw=2
    ms=8
    mew=1

    f='%spost_equal_weights.dat' % outstem
    #f='%s.txt'%outstem
    f=os.path.join(outdir,f)

    # Load equally-weighted posterior samples
    x=numpy.genfromtxt(f)#[:,2:]
    nsamp=x.shape[0]
    ncols=x.shape[1] # The fifth [seventh] column is the posterior value
    # There must be a better way, but:
    z=numpy.zeros((nsamp,ncols+nbins-1))
    z[:,:-(nbins-1)]=x
    # Shift posterior values to end
    z[:,-1]=z[:,ncols-1] # Copy...
    z[:,ncols-1]=0.0     # ...and blank

    # Fetch best-fit parameters and calculate best-fit line
    #drawmap=z[numpy.where(z[:,-1]==z[:,-1].max())]
    ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join(outdir,outstem))
    drawml=ana.get_best_fit()['parameters']

    summf=os.path.join(outdir,'1-summary.txt')
    summary=numpy.genfromtxt(summf)[-1,:]
    drawmap=summary[-(ncols+1):-2]

    #print '--> Calculating *ML* reconstruction'
    #drawmap=drawml
    power=2.5
    BM=medianArray(bins)
    NB=len(BM)

    ymap=numpy.zeros(NB)
    for ibin in xrange(NB):
        if nlaws==1:
            ymap[ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e6,drawmap[0],\
                    drawmap[1]+power,drawmap[2]/1.0e6,drawmap[3]/1.0e6,\
                    1.0/SURVEY_AREA)
        elif nlaws==2:
            ymap[ibin]=lumfunc.powerLawFuncWrap(nlaws,BM[ibin]/1.0e6,drawmap[0],\
                    drawmap[1]+power,-99.0,drawmap[4]+power,\
                    -1.0e10,1.0e10,drawmap[5]/1.0e6,-99.0,-99.0,\
                    -99.0,-99.0,1.0/SURVEY_AREA)
        elif nlaws==3:
            ymap[ibin]=lumfunc.powerLawFuncWrap(nlaws,BM[ibin]/1.0e6,drawmap[0],\
                    drawmap[1]+power,-99.0,drawmap[4]+power,\
                    -1.0e10,1.0e10,drawmap[5]/1.0e6,drawmap[6]+power,\
                    drawmap[7]/1.0e6,-99.0,-99.0,1.0/SURVEY_AREA)
        elif nlaws==4:
            ymap[ibin]=lumfunc.powerLawFuncWrap(nlaws,BM[ibin]/1.0e6,\
                    drawmap[0],drawmap[1]+power,-99.0,drawmap[4]+power,\
                    -1.0e10,1.0e10,drawmap[5]/1.0e6,drawmap[6]+power,\
                    drawmap[7]/1.0e6,drawmap[8]+power,drawmap[9]/1.0e6,\
                    1.0/SURVEY_AREA)    

    if NLAWS_SIM==0 and injectNumberOfSourcesExtracted is not None:
        NSKADS_RESCALING=373936.0/float(injectNumberOfSourcesExtracted)
        NSKADS_RESCALING=1.0
        if perGalXYRandomize and False:
            NSKADS_RESCALING=0.795

    else:
        NSKADS_RESCALING=1.0
#        NSKADS_RESCALING=0.795#1.0-0.190
        if perGalXYRandomize and False:
            NSKADS_RESCALING=0.795
    if run_num=='150121h': # HACK
        NRESCALE=0.795
    else:
        NRESCALE=1.0
        print 'Rescaling estimates for NSKADS (x %f)' % NSKADS_RESCALING
    ymap*=1.0/NRESCALE#NSKADS_RESCALING
    #print ymap

    Tb=numpy.zeros(nsamp)
    #SURVEY_AREA=2.16 # sq. deg.
    #SURVEY_AREA=1.00 # sq. deg.
    #area=REF_AREA

    # For Bondi 03, draw[0]/2.16/2. then / (4*pi) [ONLY] reconstructs ok
    # For Bondi 08, draw[0] then / (4*pi) & /S_A**2 in settings.py reconstructs ok
    #    if 'sim' in dataset:
    #    bin_medians *= 1.0e-3

    for isamp in xrange(nsamp):
        draw=z[isamp,:]
        if nlaws==1:
            Tb[isamp]=lumfunc.counts2Temperature(nlaws,draw[0],draw[1],\
                            -99.0,draw[2]/1.0e6,draw[3]/1.0e6,\
                            -99.0,-99.0,-99.0,-99.0,-99.0,\
                            freq=radioObservingFreqHz)
        elif nlaws==2:
            Tb[isamp]=lumfunc.counts2Temperature(nlaws,draw[0],draw[1],\
                                    draw[4],draw[2]/1.0e6,draw[3]/1.0e6,\
                                    draw[5]/1.0e6,-99.0,-99.0,-99.0,-99.0,\
                                    freq=radioObservingFreqHz)

        elif nlaws==3:
            Tb[isamp]=lumfunc.counts2Temperature(nlaws,draw[0],draw[1],\
                                    draw[4],draw[2]/1.0e6,draw[3]/1.0e6,\
                                    draw[5]/1.0e6,draw[6],draw[7]/1.0e6,\
                                    -99.0,-99.0,freq=radioObservingFreqHz)
        elif nlaws==4:
            Tb[isamp]=lumfunc.counts2Temperature(nlaws,draw[0],draw[1],\
                                    draw[4],draw[2]/1.0e6,draw[3]/1.0e6,\
                                    draw[5]/1.0e6,draw[6],draw[7]/1.0e6,\
                                    draw[8],draw[9]/1.0e6,\
                                    freq=radioObservingFreqHz)

            #print Tb[isamp]
        for ibin in xrange(NB):
            #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncErfsS(nlaws,1000.0*bin_medians[ibin],draw[0],draw[1],-99,-99,draw[2],draw[3],bins[ibin],bins[ibin+1],-99,SURVEY_NOISE,area)*1000.0\
            ##z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(1000.0*bin_medians[ibin],\
            ##                    draw[0],draw[1],draw[2],draw[3],area)*1000.0\
            ##                      *(1000.0*bin_medians[ibin])**2.5 \
            ##                      * sqDeg2sr / (4*pi) #/ SURVEY_AREA # ???
            S_mJy=BM[ibin] # mJy
            S_uJy=1000.0*S_mJy # mJy -> uJy
            S_Jy=S_uJy/1.0e6 # uJy -> Jy
            #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(S_Jy,\
            #                    draw[0],draw[1],draw[2]/1.0e6,draw[3]/1.0e6,SURVEY_AREA)\
            #                      *(S_Jy)**2.5 /1.4756329694112211
            #                                  * sqDeg2sr / (4*pi) #/ SURVEY_AREA # ???
            #power=2.5
            #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(S_mJy,\
            #                    draw[0]*(10**-4.5),draw[1]+power,draw[2]/1.0e3,draw[3]/1.0e3,\
            #                    SURVEY_AREA/sqDeg2sr)#*1.43982167576106

            if dataset in ['video','first','vvdf']:
                #print BM[ibin]
                # Want mJy^-1 deg^-2
                # Convert to deg^-2
                if nlaws==1:
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e6,\
                            draw[0],draw[1]+power,draw[2]/1.0e6,\
                            draw[3]/1.0e6,1.0/SURVEY_AREA)
                elif nlaws==2:
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncWrap(nlaws,\
                            BM[ibin]/1.0e6,draw[0],draw[1]+power,-99.0,draw[4]+power,\
                            #draw[2]/1.0e6,draw[3]/1.0e6,draw[5]/1.0e6,1.0/SURVEY_AREA)
                            -1.0e10,1.0e10,draw[5]/1.0e6,-99.0,-99.0,-99.0,-99.0,\
                            1.0/SURVEY_AREA)
                elif nlaws==3:
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncWrap(nlaws,\
                            BM[ibin]/1.0e6,draw[0],draw[1]+power,-99.0,draw[4]+power,\
                            #draw[2]/1.0e6,draw[3]/1.0e6,draw[5]/1.0e6,1.0/SURVEY_AREA)
                            -1.0e10,1.0e10,draw[5]/1.0e6,draw[6]+power,draw[7]/1.0e6,\
                            -99.0,-99.0,1.0/SURVEY_AREA)#/sqDeg2sr)
                elif nlaws==4:
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncWrap(nlaws,\
                            BM[ibin]/1.0e6,draw[0],draw[1]+power,-99.0,draw[4]+power,\
                            #draw[2]/1.0e6,draw[3]/1.0e6,draw[5]/1.0e6,1.0/SURVEY_AREA)
                            -1.0e10,1.0e10,draw[5]/1.0e6,draw[6]+power,draw[7]/1.0e6,\
                            draw[8]+power,draw[9]/1.0e6,1.0/SURVEY_AREA)

            elif dataset in ['cosmos']:
                # Want Jy^-1 sr^-1
                ##z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(S_uJy,\
                ##                draw[0]*(10**-1.5),draw[1]+power,draw[2],draw[3],\
                ##                SURVEY_AREA/e)#/e
                z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(\
                                BM[ibin]/1.0e3,\
                                draw[0],draw[1]+power,\
                                draw[2]/1.0e6,draw[3]/1.0e6,1.0/SURVEY_AREA)
            elif 'sim' in dataset:
                # NB BM is in uJy here
                # draw[2],draw[3] are in uJy
                # Inputs to powerLawFuncS in Jy
                #      <-> n(S) will be in Jy^1.5
                ##z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(S_uJy,\
                ##                draw[0]*(10**-1.5),draw[1]+power,draw[2],draw[3],\
                ##                SURVEY_AREA/e)#/e
                #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e3,\
                #                draw[0],draw[1]+power,draw[2]/1.0e6,draw[3]/1.0e6,\
                #                SURVEY_AREA)#/e
                #if isamp == 0 and ibin==0 and False:
                #    print draw
                #    print BM
                #    b=BM/1.0e6
                #    C=draw[0]; a=draw[1]; d2=draw[2]/1.0e6; d3=draw[3]/1.0e6
                #    print b[0],C,a,d2,d3
                #    print C*b[0]**(a+2.5)
                #    sys.exit(0)
                # Note that z samples are for SURVEY_AREA, so SURVEY_AREA -> 1.0 ?
                if nlaws==1:
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(\
                                BM[ibin]/1.0e6,\
                                draw[0],draw[1]+power,\
                                -1.0e10,1.0e10,1.0)
                                #draw[2]/1.0e6,draw[3]/1.0e6,1.0)
                #z[isamp,ncols-1+ibin]=lumfunc.C2N(\
                #                draw[0],draw[1],draw[2],draw[3],SURVEY_AREA)
                elif nlaws==2:
                    #print nlaws,S,C,alpha,D,beta,Smin,Smax,S0,area
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncDoubleS(\
                        BM[ibin]/1.0e6,\
                        draw[0],draw[1]+power,-99.0,draw[4]+power,\
                        -1.0e10,1.0e10,draw[5]/1.0e6,1.0)
                elif nlaws==3:
                    power=2.5
                    #print nlaws,S,C,alpha,D,beta,Smin,Smax,S0,area
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncWrap(nlaws,\
                        BM[ibin]/1.0e6,draw[0],draw[1]+power,-99.0,draw[4]+power,\
                        #draw[2]/1.0e6,draw[3]/1.0e6,draw[5]/1.0e6,1.0/SURVEY_AREA)
                        #draw[2]/1.0e6,draw[3]/1.0e6,draw[5]/1.0e6,draw[6]+power,draw[7]/1.0e6,\
                        -1.0e10,1.0e10,draw[5]/1.0e6,draw[6]+power,draw[7]/1.0e6,\
                        -99.0,-99.0,1.0)
                elif nlaws==4:
                    #print nlaws,S,C,alpha,D,beta,Smin,Smax,S0,area
                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncWrap(nlaws,\
                        BM[ibin]/1.0e6,draw[0],draw[1]+power,-99.0,draw[4]+power,\
                        #draw[2]/1.0e6,draw[3]/1.0e6,draw[5]/1.0e6,1.0/SURVEY_AREA)
                        -1.0e10,1.0e10,draw[5]/1.0e6,draw[6]+power,draw[7]/1.0e6,\
                        draw[8]+power,draw[9]/1.0e6,1.0)

                    #,draw[8]+power,draw[9]/1.0e6,1.0)
            #elif 'first' in dataset:
                # NOW DEALT WITH ABOVE
                # units are:   draw[0] sr^-1 Jy^-1
                #                  [1] none
                #                  [2] uJy
                #                  [3] uJy
                #
                #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(S_uJy,\
                #        draw[0],draw[1]+power,draw[2],draw[3],\
                #        SURVEY_AREA/sqDeg2sr)#/e
                # HACK SMIN/SMAX FOR PLOTTING - REDO THIS
            #    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e6,\
            #                    draw[0],draw[1]+power,\
            #                    draw[2]/1.0e6,draw[3]/1.0e6,SURVEY_AREA)
                #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e6,\
                #                7.2,-2.11+power,\
                #                draw[2]/1.0e6,draw[3]/1.0e6,SURVEY_AREA)
                #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e6,\
                #                19.7/SURVEY_AREA,-2.32+power,\
                #                draw[2]/1.0e6,draw[3]/1.0e6,SURVEY_AREA)

                
            #z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncS(S_mJy,\
                                #draw[0],draw[1],draw[2]/1.0e3,draw[3]/1.0e3,\
                                #SURVEY_AREA/sqDeg2sr)*(S_mJy)**2.5

            #print S_Jy, draw[0], draw[1], draw[2], draw[3]
                                #SURVEY_AREA, lumfunc.powerLawFuncS(S_Jy,draw[0],draw[1],\
                                #              draw[2]/1.0e6,draw[3]/1.0e6,SURVEY_AREA)\
                                #              *S_Jy**2.5,\
                                #              z[isamp,ncols-1+ibin]
                                #raw_input('Press enter to continue')

    if False and perGalXYRandomize:
        NSKADS_RESCALING=1.0#0.795
        z*=1.0/NSKADS_RESCALING#/0.795#NSKADS_RESCALING#*sqDeg2sr
    if run_num=='150121h': # HACK
        z*=1.0/NRESCALE

    #z*=373936.0/float(injectNumberOfSourcesExtracted)
    #z*=72000.0/374995.0

    z[numpy.where(z==0.0)]='NaN' # 0.0 -> NaN blanking
    #print BM*1000.0

    # Save the raw reconstructions
    reconf='recon_raw.txt'
    reconf=os.path.join(outdir,reconf)
    recons=z[:,ncols-1:-1]
    numpy.savetxt(reconf,recons)

    # Save the brightness temperatures
    tempsf='temperatures.txt'
    tempsf=os.path.join(outdir,tempsf)
    Tb = Tb[~numpy.isnan(Tb)]
    numpy.savetxt(tempsf,Tb)
    #tt=calculate_confidence(Tb,alpha=0.68)
    print 'Tb / mK = %f +/- %f' % (1.0e3*numpy.mean(Tb),1.0e3*numpy.std(Tb)/numpy.sqrt(len(Tb)))

    # Generate stats here...
    s=numpy.zeros((nbins-1,6))
    s[:,0]=BM
    # Fix this sim special case later:
    #if 'sim' not in dataset: s[:,0]*=1000.0 # mJy -> uJy
    print '# ibin flux fit low high dlower dupper skew kurtosis'
    for ibin in xrange(nbins-1):
        x = recons[:,ibin]
        # Remove NaNs from stats vectors
        # http://stackoverflow.com/questions/11620914/removing-nan-values-from-an-array
        x = x[~numpy.isnan(x)]
        #ss=stats.bayes_mvs(x,alpha=0.68)[0]
        #x*=numpy.power(s[ibin,0]/1.0e6,2.5)

        try:
            ss=numpy.zeros(3)
            ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence(x,alpha=0.68,ret_all=True)
        except:
            ss=numpy.nan*numpy.ones(3)
        tt=peak_confidence(x,bins=10)
        #tt*=72000.0/374995.0
        #ss*=72000.0/374995.0
        #ss*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print ss[0],tt
        s[ibin,1]=ss[0]  # median
        #s[ibin,1]=tt     # peak
        #s[ibin,1]=ymap[ibin] # MAP
        #print ymap
        #sys.exit(0)
        s[ibin,2]=ss[1]  # lower
        s[ibin,3]=ss[2]  # upper
        s[ibin,4]=stats.skew(x) # skewness
        s[ibin,5]=stats.kurtosis(x) # kurtosis
     #   if perGalXYRandomizeMask:
     #       s[ibin,1:]*=1.0/perGalXYRandomizeMaskArea
        print ibin,s[ibin,0],s[ibin,1],dlow,dhigh,ss[1],ss[2],s[ibin,4],s[ibin,5]#,stats.skewtest(x)

    # ...and output to file
    rstatsf='recon_stats.txt'
    rstatsf=os.path.join(outdir,rstatsf)
    hdr='# median_flux_uJy dnds_2p5_Jy1p5srm1 delta_dnds_2p5_lower_Jy1p5srm1 delta_dnds_2p5_upper_Jy1p5srm1 skewness kurtosis'
    fid = open(rstatsf,'w')
    print hdr
    print s

    print '-> Writing stats (i.e. dnds2p5+/- v. S/uJy) to %s' % rstatsf
    #print '(you can plot these in topcat alongside 1000.0*corr*bondi200X_orig.txt)'
    fid.write('%s\n'%hdr)
    numpy.savetxt(fid,s)
    fid.close()

#-------------------------------------------------------------------------------
# Now do the plotting:
#=====================

    # Fetch and plot input datafile
    if dataset in ['video']:
        datafile='video/%s_test_%i_%s.txt'%(perGalGalaxyType,nb,run_num_run)
    elif dataset in ['sim']: datafile=os.path.join(dataset,datafile)
    #datafile='sims/141003d/sim.txt'
    #datafile='sims/141013a/sim_extracted.txt'
    #datafile='sims/141009c/sim.txt'
    print '-> Plotting: datafile = %s' % datafile
    true=numpy.genfromtxt(datafile)
    # Strip out any negative bins
    true=true[numpy.where(true[:,0]>=0.0)]
    #if 'sim' in datafile:
    xb=true[:,2]      # uJy
    #else:
    #    xb=1000.0*true[:,2]      # mJy -> uJy
    #if 'video' in datafile:
    #    xb=true[:,2]
    # C factors:
    # 1.05 - resolution - done by binner.py
    # 0.88 - completeness - done by CORR_BINS in binner.py
    # 0.97 Halo area factor - done here (usually)
    S_A=1.00 # This is not needed - the binfile is already in standard units
    #S_A=0.795
    #if perGalXYRandomizeMask:
    #    S_A=perGalXYRandomizeMaskArea
    print '***Warning - setting SURVEY_AREA = %4.2f TEMPORARILY' % S_A
#    true[:,5:7]*=373936.0/float(injectNumberOfSourcesExtracted)
    if dataset in ['cosmos','vvdf']:
        xb=1000.0*true[:,2] # mJy -> uJy
    yb=true[:,5]*true[:,8] / S_A  #/ SURVEY_AREA # Apply the correction factor
    dyb=true[:,6]*true[:,8] / S_A #/ SURVEY_AREA # Apply the correction factor
    print '-> Corrected counts for completeness and/or ang. size (C factor)'
    if dyb[0]==-99: dyb=0.0  # handle sims
    #print yb
    #sys.exit(0)
    #if 'sim' in datafile: # ?
    #    # In the sim files, dn/ds is in standard units
    #    # Convert to SURVEY_AREA^-1
    #    # Convert to Jy^1.5, i.e. 
    #    yb *= sqDeg2sr*SURVEY_AREA
    #    dyb *= sqDeg2sr*SURVEY_AREA
    #print xb
    #print yb
    #yb*=1.0/(1000.0*sqDeg2sr) ############# JUSTIFY ###########
    #yb*=1.0*0.876#/5.25 # Simulation area hack
    #print 'Hacked simulation area for 8x to 1/0.876 [again]'
    df='k+-.' #'b+-.'
    print xb
    print yb
    #sys.exit(0)
    #plt.errorbar(xb[:-1],yb[:-1],yerr=dyb[:-1],fmt=df,capsize=0,label='data (extracted)')
    if 'sim' not in dataset and (S0_MIN!=S0_MAX) and (SMIN_SKADS==0.01):
        plt.errorbar(xb,yb,yerr=dyb,fmt=df,capsize=0,label='data (extracted)',\
                     lw=lw,ms=ms,mew=mew)

    if 'sim' in datafile: # Hack to plot noise-free SKADS data
        #true_noise_free=numpy.genfromtxt('sims/141003d/sim.txt')
        noisefreef='sims/%s/sim_noiseless.txt' % run_num_run
        true_noise_free=numpy.genfromtxt(noisefreef)
        true_noise_free=true_noise_free[numpy.where(true_noise_free[:,0]>=0.0)]
        xb=true_noise_free[:,2]
        NSKADS_RESCALING=1.0
        true_noise_free[:,5:7]*=NSKADS_RESCALING
        yb=true_noise_free[:,5]*true_noise_free[:,8]
        dyb=true_noise_free[:,6]*true_noise_free[:,8]
        if dyb[0]==-99: dyb=0.0
        df='g+--'
        plt.errorbar(xb,yb,yerr=dyb,fmt=df,capsize=0,label='data (noise-free, input)',lw=lw,ms=ms,mew=mew)

    # For FIRST/VIDEO/VVDF, also plot the COSMOS data....
    if dataset in ['first','video','vvdf']:
        if dataset != 'vvdf':
            cosdataf='cosmos/bondi2008_%s.txt' % 'all'
        else:
            cosdataf='vvdf/bondi2003_orig.txt'

        print '-> Plotting: datafile = %s' % cosdataf
        cos=numpy.genfromtxt(cosdataf)
        xc=1.0e3*cos[:,2]      # mJy -> uJy
        ycc=cos[:,5]*1.0  #/ SURVEY_AREA # Don't apply the correction factor
        dycc=cos[:,6]*1.0  #/ SURVEY_AREA # Don't apply the correction factor
        if False:
            plt.errorbar(xc,ycc,yerr=dycc,fmt='g*',fillstyle='none',markerfacecolor='w',\
                     capsize=0,\
                     label='%s uncorrected (%s)'%(cosdataf.split('/')[0].upper(),\
                                                  cosdataf.split('/')[-1].split('_')[0]))

        yc=cos[:,5]*cos[:,8]  #/ SURVEY_AREA # Apply the correction factor
        dyc=cos[:,6]*cos[:,8]  #/ SURVEY_AREA # Apply the correction factor
        print '-> Corrected counts for completeness and/or ang. size (C factor)'
        if False:
            plt.errorbar(xc,yc,yerr=dyc,fmt='g*',markeredgecolor='g',capsize=0,\
                     label='%s corrected (%s)'%(cosdataf.split('/')[0].upper(),\
                                                cosdataf.split('/')[-1].split('_')[0]))

    # ... and the VVDF data
    if dataset in ['video']:
        cosdataf='vvdf/bondi2003_%s.txt' % 'all'
        print '-> Plotting: datafile = %s' % cosdataf
        cos=numpy.genfromtxt(cosdataf)
        xc=1.0e3*cos[:,2]      # mJy -> uJy
        ycc=cos[:,5]*1.0  #/ SURVEY_AREA # Don't apply the correction factor
        dycc=cos[:,6]*1.0  #/ SURVEY_AREA # Don't apply the correction factor

        old_label=cosdataf.split('/')[-1].split('_')[0]
        plt.errorbar(xc,ycc,yerr=dycc,fmt='b*',fillstyle='none',markerfacecolor='w',
                     capsize=0,\
                     label='%s uncorrected (%s)'%(cosdataf.split('/')[0].upper(),\
                                                  'Bondi+ 2003'))

        yc=cos[:,5]*cos[:,8]  #/ SURVEY_AREA # Apply the correction factor
        dyc=cos[:,6]*cos[:,8]  #/ SURVEY_AREA # Apply the correction factor
        print '-> Corrected counts for completeness and/or ang. size (C factor)'
        old_label=cosdataf.split('/')[-1].split('_')[0]
        plt.errorbar(xc,yc,yerr=dyc,fmt='b*',markeredgecolor='b',capsize=0,\
                     label='%s corrected (%s)'%(cosdataf.split('/')[0].upper(),\
                                                'Bondi+ 2003'))


    # and Kim's data.....
    if dataset in ['video']:
        cosdataf='mca/bondi2003_mca.txt'
        print '-> Plotting: datafile = %s' % cosdataf
        cos=numpy.genfromtxt(cosdataf)
        xc=cos[:,2]      # mJy -> uJy
        S_AREA=1.0
        ycc=cos[:,5]*S_AREA  #/ SURVEY_AREA # Don't apply the correction factor
        dycc=cos[:,6]*S_AREA  #/ SURVEY_AREA # Don't apply the correction factor
        if False:
            plt.errorbar(xc,ycc,yerr=dycc,fmt='r*',fillstyle='none',markerfacecolor='w',
                     capsize=0,\
                     label='%s uncorrected (%s)'%(cosdataf.split('/')[0].upper(),\
                                                  cosdataf.split('/')[-1].split('_')[0]))

        yc=cos[:,5]*cos[:,8]*S_AREA  #/ SURVEY_AREA # Apply the correction factor
        dyc=cos[:,6]*cos[:,8]*S_AREA  #/ SURVEY_AREA # Apply the correction factor
        print '-> Corrected counts for completeness and/or ang. size (C factor)'
        if False:
            plt.errorbar(xc,yc,yerr=dyc,fmt='r*',markeredgecolor='r',capsize=0,\
                     label='%s corrected (%s)'%(cosdataf.split('/')[0].upper(),\
                                                cosdataf.split('/')[-1].split('_')[0]))

    # Plot the reconstruction
    offset=POINTS_OFFSET # 0.0 was 5.0
    xrecon=s[:,0]
    if 'cosmos' in datafile:
        xrecon=1000.0*s[:,0]
    yrecon=s[:,1]
    #if 'first' in datafile: yrecon *= 1.0/SURVEY_AREA # ???? NEED TO JUSTIFY THIS
#    s[:,1:]*=1.0*0.876#/5.25 # Simulation area hack
    #s[:,1:]*=37.5
#    print 'Hacked simulation area for 8x to 1/0.876'
    #yrecon_down=s[:,2]; yrecon_up=s[:,3]
    yrecon_down=s[:,2]-s[:,1]; yrecon_up=s[:,3]-s[:,1]
    if False:
        plt.errorbar(xrecon+offset,yrecon,fmt='c.',alpha=1.0,label='median reconstruction',\
                 lw=lw,ms=ms,mew=mew)
    #plt.fill_between(xrecon+offset,yrecon-dyrecon_down,yrecon+dyrecon_up,\
    #                 color='k',alpha=0.2)
    if True:
        plt.fill_between(xrecon+offset,ymap+yrecon_down,ymap+yrecon_up,\
                     color='k',alpha=0.2)
        plt.errorbar(xrecon+offset,ymap,fmt='k',label='MAP estimate',\
                     lw=lw,ms=ms,mew=mew)

    #recons*=numpy.power(xrecon/1.0e6,2.5)
    #xrecons=numpy.tile(xrecon+offset,len(recons[:,0]))
    #plt.scatter(xrecons,recons,alpha=0.01,c='k',s=2,edgecolors='None')

    #plt.errorbar(xrecon+offset,yrecon-dyrecon_down,fmt='k-',alpha=0.2)
    #plt.errorbar(xrecon+offset,yrecon+dyrecon_up,fmt='k-',alpha=0.2)
    #plt.errorbar(xrecon+offset,yrecon,yerr=[dyrecon_down,dyrecon_up],\
    #             fmt='k+-',capsize=0,label='reconstruction')
    # Do a quick fit to the reconstructed points
    #print xrecon
    #print yrecon
    xr=numpy.log10(xrecon/1.0e6); yr=numpy.log10(yrecon)
    xr = xr[numpy.isfinite(xr)]; yr = yr[numpy.isfinite(yr)]
    g=numpy.polyfit(xr,yr,1)
    print '-> Fitted power law =~ %4.2f . S ^%4.2f' % (10.0**g[1],g[0]-power)
    #for isamp in xrange(nsamp):
    #    plt.errorbar(xrecon+offset,recons[isamp,:],fmt='k',alpha=0.001)
    #plt.errorbar(xrecon+offset,yrecon,fmt='k',alpha=1.0,label='reconstruction')

    # Calculate confusion noise for this model
    xinterp=numpy.linspace(xrecon[0],xrecon[-1],1001)
    yinterp=numpy.interp(xinterp,xrecon,yrecon)
    sigmaConf_jz=numpy.sqrt(lumfunc.calculateConfusionNoiseSqArray(radioSynthOmegaSr,\
            xinterp/1.0e6,yinterp*numpy.power(xinterp/1.0e6,-2.5),-99.0,5.0*SURVEY_NOISE/1.0e6))
    print 'sigmaConf_jz / uJy = %f' % (sigmaConf_jz*1.0e6)

    # Calculate confusion noise as fn of SMAX
    if True:
        smaxes=numpy.linspace(0.1,5.0,100)
        conff='confusion-noise.txt'
        conff=os.path.join(outdir,conff)
        confs=numpy.zeros((len(smaxes),2))
        #print '# smax confusion_noise'
        iS=0
        for smax_frac in smaxes:
            confs[iS,0]=smaxes[iS]
            confs[iS,1]=1.0e6*numpy.sqrt(lumfunc.calculateConfusionNoiseSqArray(radioSynthOmegaSr,xinterp/1.0e6,yinterp*numpy.power(xinterp/1.0e6,-2.5),-99.0,smax_frac*SURVEY_NOISE/1.0e6))/SURVEY_NOISE
            iS+=1
        numpy.savetxt(conff,confs)
        print '--> Post-process %s with confusion-noise.ipynb' % conff
    

    # Sanity check via an array
    fl,co=lumfunc.powerLaw2Array(1,10.0**g[1],g[0]-power,-99.0,-99.0,-99.0,-99.0,-99.0,-99.0,xinterp[0]/1.0e6,xinterp[-1]/1.0e6)
    sigmaConf_jz2=numpy.sqrt(lumfunc.calculateConfusionNoiseSqArray(radioSynthOmegaSr,\
            fl,co,-99.0,5.0*SURVEY_NOISE/1.0e6))
    print 'sigmaConf_jz2 / uJy = %f' % (sigmaConf_jz2*1.0e6)

    #if nlaws==1:
    #    sigmaConf_jz=numpy.sqrt(lumfunc.calculateConfusionNoiseSq(radioSynthOmegaSr,nlaws,10.0**g[1],g[0]-power,-99.0,-99.0,-99.0,-99.0,50.0e-9,5.0*SURVEY_NOISE/1.0e6))
    #    print 'sigmaConf_jz / uJy = %f' % (sigmaConf_jz*1.0e6)

    # Calculate Condon-style confusion figure for this model
    # SLOW!
    if False:
        sigmaConf_star_jz=pofd.dnds2conf(xinterp/1.0e6,yinterp*numpy.power(xinterp/1.0e6,-2.5),\
                                         radioSynthBeamFWHM*radioPixels2Arcsec,radioPixels2Arcsec)
        print 'sigmaConf_star_jz / uJy = %f' % sigmaConf_star_jz

    # Plot the original power law (if it exists!)
    if 'video' in dataset:
        xorig=BM
        yorig=numpy.zeros(len(xorig))
        yorig2=numpy.zeros(len(xorig))
        yorig3=numpy.zeros(len(xorig))
        dyorig=-99
        for ibin in xrange(nbins-1):
            #C_TRUE=19.7; ALPHA_TRUE=-2.32
            C_TRUE=57.54/2.107; ALPHA_TRUE=-2.28; AREA_TRUE=0.97
            yorig[ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e6,\
                                              C_TRUE,\
                                              ALPHA_TRUE+power,\
                                              SMIN_TRUE/1.0e6,\
                                              SMAX_TRUE/1.0e6,AREA_TRUE)
            yorig[ibin]=lumfunc.powerLawFuncWrap(2,BM[ibin]/1.0e6,0.691126368989902112E+05,-1.50+power,-99.0,-1.74+power,0.237/1.0e6,265.3/1.0e6,58.5/1.0e6,-99.0,-99.0,-99.0,-99.0,0.97)
            yorig[ibin]=lumfunc.powerLawFuncWrap(3,BM[ibin]/1.0e6,0.127E+04,-1.72+power,-99.0,-0.15+power,0.80e-1/1.0e6,72.4/1.0e6,4.4/1.0e6,-2.47+power,20.5/1.0e6,-99.0,-99.0,0.97)
        #plt.errorbar(xorig,yorig,yerr=dyorig,fmt='b--',capsize=0,\
        #             label='$%6.2fS^{%6.2f} . S^{2.5}$ (bondi2003)'%(C_TRUE,ALPHA_TRUE))

    elif 'video' not in dataset:
        xorig=BM # uJy
        if dataset in ['cosmos']:
            xorig=1000.0*BM # mJy -> uJy
        yorig=numpy.zeros(len(xorig))
        yorig2=numpy.zeros(len(xorig))
        yorig3=numpy.zeros(len(xorig))
        dyorig=-99
        #print xorig,C_TRUE,ALPHA_TRUE,power,SMIN_TRUE,SMAX_TRUE,SURVEY_AREA

        for ibin in xrange(nbins-1):
            #print BM[ibin],bins[ibin],bins[ibin+1]
            if dataset in ['first','sim','vvdf'] or 'sim' in dataset:
                if dataset in 'first':
                    C_TRUE=19.7; ALPHA_TRUE=-2.32
                elif 'sim' in dataset:
                    C_TRUE=C_SIM; ALPHA_TRUE=ALPHA_SIM
                    CONVERT_C_TRUE=10.0**(6.0*(ALPHA_TRUE+2.5))
                    print C_TRUE,CONVERT_C_TRUE,ALPHA_TRUE
                    C_TRUE=CONVERT_C_TRUE/C_TRUE
                    BETA_TRUE=BETA_SIM
                    S0_TRUE=S0_SIM
                elif dataset in 'vvdf':
                    C_TRUE=57.54/2.107; ALPHA_TRUE=-2.28
                yorig[ibin]=lumfunc.powerLawFuncS(xorig[ibin]/1.0e6,\
                                                  C_TRUE,\
                                                  ALPHA_TRUE+power,\
                                                  SMIN_TRUE/1.0e6,\
                                                  SMAX_TRUE/1.0e6,1.0)

                if NLAWS_SIM==2:
                    yorig[ibin]=lumfunc.powerLawFuncDoubleS(xorig[ibin]/1.0e6,\
                                                    C_TRUE,ALPHA_TRUE+power,-99.0,\
                                                    BETA_TRUE+power,\
                                                    SMIN_TRUE/1.0e6,\
                                                    SMAX_TRUE/1.0e6,
                                                    S0_TRUE/1.0e6,1.0)

#                    z[isamp,ncols-1+ibin]=lumfunc.powerLawFuncDoubleS(\
#                        BM[ibin]/1.0e6,\
#                        draw[0],draw[1]+power,-99.0,draw[4]+power,\
#                        draw[2]/1.0e6,draw[3]/1.0e6,draw[5]/1.0e6,1.0)
                C_TRUE2=7.2; ALPHA_TRUE2=-2.11
                #C_TRUE2=C_SIM/sqDeg2sr; ALPHA_TRUE2=ALPHA_SIM
                yorig3[ibin]=lumfunc.powerLawFuncS(BM[ibin]/1.0e6,\
                                          C_TRUE2,\
                                          ALPHA_TRUE2+power,\
                                          SMIN_TRUE/1.0e6,SMAX_TRUE/1.0e6,1.0)
            else:
                # NEED TO CHECK FACTORS OF SURVEY_AREA HERE
                yorig[ibin]=lumfunc.powerLawFuncS(xorig[ibin]/1.0e6,\
                                          C_TRUE/SURVEY_AREA,ALPHA_TRUE+power,\
                                          SMIN_TRUE/1.0e6,SMAX_TRUE/1.0e6,SURVEY_AREA)

#            yorig2[ibin]=lumfunc.powerLawFuncErfsS(BM[ibin]/1.0e6,\
#                    NLAWS_SIM,C_TRUE,ALPHA_TRUE+power,D_SIM,BETA_SIM,\
#                    SMIN_TRUE/1.0e6,SMAX_TRUE/1.0e6,\
#                    SMIN_TRUE/1.0e6,SMAX_TRUE/1.0e6,\
#                    S0_SIM/1.0e6,GAMMA_SIM+power,S1_SIM/1.0e6,\
#                    DELTA_SIM+power,S2_SIM/1.0e6,NOISE_SIM/1.0e6,SURVEY_AREA)
                    
        #yorig *= 1.0/sqDeg2sr
        #print yorig*sqDeg2sr
        print yb
        print yorig
        print yrecon
        #print yb/yorig,yrecon/yorig#,yb/yorig2,yrecon/yorig2,yorig/yorig2
#        plt.errorbar(xorig,yorig,yerr=dyorig,fmt='g--',capsize=0,\
#                     label='$%6.2fS^{%6.2f} . S^{2.5}$'%(C_TRUE,ALPHA_TRUE))
        #plt.errorbar(xorig,yorig3,yerr=dyorig,fmt='m--',capsize=0,\
        #             label='$%6.2fS^{%6.2f} . S^{2.5}$ (KMW14)'%(C_TRUE2,ALPHA_TRUE2))

    # Overlay literature values
    litfile='dezotti/1d4GHz.dat'
    lit=numpy.genfromtxt(litfile)
    xlit=lit[:,0]; ylit=lit[:,1]
    dylit_up=lit[:,2]; dylit_down=lit[:,3]
    if 'sim' not in dataset:
        plt.errorbar(10**(6+xlit),ylit,yerr=[dylit_up,dylit_down],\
                 fmt='k.',alpha=0.25,label='literature')

    # Overlay SKADS values
    skadsfile='heywood/wilman_counts_1p4GHz_jz.txt'
    skads=numpy.genfromtxt(skadsfile)
    xskads=1.0e6*skads[:,0]
    yskads=skads[:,3]

    #NSKADS_RESCALING=71962.0/373936.0
    #yskads *= NSKADS_RESCALING

    #dyskads=skads[:,4]
    if True:
        alpha=0.5 # 0.25
        plt.errorbar(xskads,yskads,fmt='b',alpha=alpha,label='Wilman+ 2008',\
                                                       lw=lw,ms=ms,mew=mew)

    # Overlay Vernstrom values
    vernfile='tessa/vernstrom_model_1p4_jz.dat'
    vern=numpy.genfromtxt(vernfile)
    xvern=vern[:,0]
    yvern_down= 10**(2.5*xvern + vern[:,1])
    yvern=      10**(2.5*xvern + vern[:,2])
    yvern_up=   10**(2.5*xvern + vern[:,3])
    if 'sim' not in dataset:
        alpha=0.5 # 0.25
        plt.errorbar(10**(6+xvern),yvern,fmt='r',alpha=alpha,label='Vernstrom+ 2014',\
                                                             lw=lw,ms=ms,mew=mew)
        plt.errorbar(10**(6+xvern),yvern_down,fmt='r-.',alpha=alpha,lw=lw,ms=ms,mew=mew)
        plt.errorbar(10**(6+xvern),yvern_up,fmt='r-.',alpha=alpha,lw=lw,ms=ms,mew=mew)
    # Calculate temperature contribution for Vernstrom values above some Smin
    Tb_vernstrom=lumfunc.counts2TemperatureArray(radioSynthOmegaSr,
                    10**vern[:,0],10**vern[:,2],-99.0,5.0*SURVEY_NOISE/1.0e6,\
                    freq=radioObservingFreqHz)
    print 'Tb_vernstrom / mK = %f' % (Tb_vernstrom*1.0e3)

    # Calculate confusion noise contribution for Vernstrom values, Smin<S<Smax
    if False:
        sigmaConf_v=numpy.sqrt(lumfunc.calculateConfusionNoiseSqArray(radioSynthOmegaSr,\
                    10**vern[:,0],10**vern[:,2],-99.0,5.0*SURVEY_NOISE/1.0e6))
        print 'sigmaConf_vernstrom / uJy = %f' % (sigmaConf_v*1.0e6)
        sigmaConf_star_v=pofd.dnds2conf(10.0**vern[:,0],10.0**vern[:,2],radioSynthBeamFWHM*radioPixels2Arcsec,radioPixels2Arcsec)
        print 'sigmaConf_vernstrom / uJy = %f' % sigmaConf_star_v
        #sys.exit(0)

    # Take care of the furniture and plot to file
    f=24
    plotf='recon_plot_%s.pdf' % run_num
    if False:
        plt.title('%s - %s'%(outdir,dataset))
    plt.xscale('log')
    plt.yscale('log')
    #if 'video' not in dataset:
    plt.xlim(10.0,1000.0) # xx lim
    if 'sim' in dataset or 'video' in dataset:
            plt.xlim(1.0,10000.0) # xx lim
            #plt.xlim(50.0,10000.0)

    #if 'first' in dataset:
    #    plt.ylim(0.06,50.0) # yy lim
    #elif 'video' in dataset:
    plt.ylim(0.06,50.0) # yy lim
    #else:
    #    pass
    #plt.ylim(0.1,100.0) # yy lim
    if 'video' in dataset:
        plt.ylim(0.05,50.0)
        #plt.ylim(1.0,50.0)
    plt.xlim(0.01,1.0e4)
    plt.ylim(1.0e-4,50.0)

    #plt.xlim(-50.0,50.0)

    # Show 1- and 5-sigma limits
    plt.axvline(SURVEY_NOISE,color='b',alpha=0.2,linewidth=2)
    plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2,linewidth=2)
    plt.text(17.0,0.05e-2,'1 sigma',rotation=90,color='b',alpha=0.5,size=12)
    plt.text(85.0,0.05e-2,'5 sigma',rotation=90,color='b',alpha=0.5,size=12)
    if 'sim' in dataset and (S0_MIN==S0_MAX):
        #ARROW_WIDTH={10.0:0.15,5.0:0.08,1.0:0.015,0.5:0.008}
        #width=ARROW_WIDTH[S0_MIN]; head_width=width*10.0
        #plt.arrow(S0_MIN,1.0e-4,0.0,3.0e-4,alpha=1.0,color='r',width=width,ec='r',\
        #          head_width=head_width,head_length=1.0e-4)
        plt.axvline(S0_MIN,color='r',ls='dashed',alpha=1.0,linewidth=2)
    if 'sim' in dataset and (SMIN_SKADS > 0.01):
        #ARROW_WIDTH={0.1:0.0015,0.5:0.01,1.0:0.02}
        #width=ARROW_WIDTH[SMIN_SKADS]; head_width=width*10.0
        #plt.arrow(SMIN_SKADS,1.0e-4,0.0,3.0e-4,alpha=1.0,color='r',width=width,ec='r',\
        #          head_width=head_width,head_length=1.0e-4)
        #plt.axvline(SMIN_SKADS,color='r',alpha=0.2,linewidth=2)
        pass
    plt.text(7.0e3,2.0e-4,PLOT_LABEL,color='k',size=48,horizontalalignment='right')

    plt.legend(loc='upper left',prop={'size':12},frameon=False,numpoints=1)
    plt.xlabel('$S/\mu\mathrm{Jy}$',fontsize=f)
    plt.ylabel('$nS^{2.5}/\mathrm{sr}^{-1}\mathrm{Jy}^{1.5}$',fontsize=f)

    plotf=os.path.join(outdir,plotf)
    plt.savefig(plotf)
    plt.close()
    print '-> Run: open %s' % plotf


    if 'sim' in dataset and False:
        fluxes=numpy.genfromtxt(os.path.join(dataset,'fluxes.txt'))
        counts=numpy.histogram(fluxes,bins=bins_saved)[0]
        plt.errorbar(medianArray(bins_saved),counts,yerr=numpy.sqrt(counts),fmt='r-.',label='simulation (noiseless)')
        print counts.sum()
        fluxes_noisy=numpy.genfromtxt(os.path.join(dataset,'fluxes_noisy.txt'))
        counts=numpy.histogram(fluxes_noisy,bins=bins_saved)[0]
        plt.errorbar(medianArray(bins_saved),counts,yerr=numpy.sqrt(counts),fmt='k.',label='simulation')

        noise=counts.sum()*gaussian(medianArray(bins_saved),0.0,SURVEY_NOISE,norm=True)
        plt.errorbar(medianArray(bins_saved),noise,fmt='g',label='noise')
        print counts.sum()
        recon_law=7236.34*numpy.power(medianArray(bins_saved),-1.72+1)
        plt.errorbar(medianArray(bins_saved),recon_law,fmt='b',label='fit')

        bins_tail=medianArray(bins_saved)[numpy.where(medianArray(bins_saved)<=0.0)]
        bins_tail=numpy.concatenate((bins_tail,abs(bins_tail[::-1])))
        counts_noise=medianArray(counts)[numpy.where(medianArray(bins_saved)<=0.0)]
        counts_noise=numpy.concatenate((counts_noise,counts_noise[::-1]))
        print bins_tail
        print counts_noise
        
        plt.xlim(bins_saved[0],bins_saved[-1])
        plt.ylim(1.0e2,1.0e5)

        f=24
        plt.title('%s - %s'%(outdir,dataset))
        plt.yscale('log')

        plt.legend(loc='upper right',prop={'size':8},frameon=False,numpoints=1)
        plt.xlabel('$S/\mu\mathrm{Jy}$',fontsize=f)
        plt.ylabel('counts')

        histof='recon_histo.pdf'
        histof=os.path.join(outdir,histof)
        plt.savefig(histof)
        plt.close()
        print '-> Run: open %s' % histof
    
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
