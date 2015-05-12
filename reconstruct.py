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
setf='%s.bayestack_settings' % param_file

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

    for isamp in xrange(nsamp):
        draw=z[isamp,:]
        for ibin in xrange(NB):
            S_mJy=BM[ibin] # mJy
            S_uJy=1000.0*S_mJy # mJy -> uJy
            S_Jy=S_uJy/1.0e6 # uJy -> Jy
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

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
