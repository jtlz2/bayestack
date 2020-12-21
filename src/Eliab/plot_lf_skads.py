#!/usr/bin/env python

"""
This is plot_counts.py
Jonathan Zwart
May 2014

Usage:

./plot_counts.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
if os.getenv('PBS_O_HOST') not in ['baltasar']:
    from matplotlib import pyplot as plt
import pylab
from profile_support import profile
import pymultinest
import lumfuncUtils
from bayestackClasses import countModel
from utils import sqDeg2sr,fetchStats
from countUtils import calculateDnByDs,medianArray,powerLawFuncQuadS
from lumfuncUtils import get_Lbins,get_sbins,get_dl, get_Vmax,get_dsdl,get_z
from plat import open_anytype, calc_Vmax, calcu_zeff,lumfuncs
from matplotlib.ticker import AutoMinorLocator
from cycler import cycler
#import sourcecounts1

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file

@profile
def main():

    """
    """

    print 'Settings file is %s' % setf

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
    recon_expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
    #print expt.data
    #print expt.bins

    # Get MAP parameters
    print '\n'*5
    print len(expt.parameters)
    ncols = len(expt.parameters)
    summf=os.path.join(outdir,'1-summary.txt')
    summary=numpy.genfromtxt(summf)[-1,:]
    drawmap=summary[-(ncols+2):-2] 
    ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join(outdir,outstem))
    drawmap=ana.get_best_fit()['parameters']
    
    if dataset=='sdss':
        print datafiles, len(datafiles[0])
        dataf = datafiles[0]
        if len(dataf)==20:
            num = dataf[-5]
            print num
        else:
            num = dataf[-6:-4]
        d=numpy.genfromtxt('%s/%s'%(outdir,datafiles[0][4:]))
        print datafiles[0][4:]
    elif dataset=='cosmos':
        print datafiles, len(datafiles[0])
        dataf = datafiles[0]
        if len(dataf)==24:
            num = dataf[-5]
            print num
        else:
            num = dataf[-6:-4]
        d=numpy.genfromtxt('%s/%s'%(outdir,datafiles[0][8:]))
        print datafiles[0][8:]

    elif 'sim' in dataset:
        d=numpy.genfromtxt(os.path.join(outdir,'sim.txt'))
        
    ssbin = d[:,2]
    Nbin  = d[:,3]
    dnds = d[:,4]    
    sbin2 = ssbin
    dnds2 = dnds 
    #bins2 = numpy.arange(21.4,29.2,0.4)
    bins2 = numpy.arange(18.0,29.2,0.4)
    bins3 = numpy.arange(18.,29.2,0.01)
    print len(expt.parameters)
    print 'this is 1st num ', num
    #sys.exit()
    params = drawmap# [0]
    #for i in drawmap:
    #	params.append(i)
    print params
    print expt.parameters

    if modelFamily in ['LFsch','LFdpl_pl','LFdpl_dpl','LFlognorm_dpl']:
        Lnorm=params[expt.parameters.index('LNORM')]
        Lstar=params[expt.parameters.index('LSTAR')]
        Lslope=params[expt.parameters.index('LSLOPE')]
        #Lzevol=params[expt.parameters.index('LZEVOL')]
    if modelFamily in ['LFdpl_pl','LFdpl_dpl','LFlognorm_dpl']:
        Lslope2=params[expt.parameters.index('LSLOPE2')]
    if modelFamily in ['LFlognorm','LFpl', 'LFdpl', 'LFdpl_pl', 'LFlognorm_dpl','LFdpl_dpl']:
        Lnorm_2=params[expt.parameters.index('LNORM_2')]
        Lstar_2=params[expt.parameters.index('LSTAR_2')]
        Lslope_2=params[expt.parameters.index('LSLOPE_2')]

    if modelFamily in ['LFdpl_dpl','LFdpl']:
        Lslope2_2=params[expt.parameters.index('LSLOPE2_2')]
    
    if modelFamily in ['LFlognorm','LFlognorm_dpl']:
        Lsigma = params[expt.parameters.index('LSIGMA')]

    Lmin=params[expt.parameters.index('LMIN')]
    Lmax=params[expt.parameters.index('LMAX')]
    truth= {'noise': 150,'LMIN':22.,'LMAX':24.5,'LNORM':numpy.log10(1e-7),'LSTAR': 23.3,'LSLOPE':3.,'LSLOPE2':1.}    
    #z = [0, 0.5, 1 , 1.5, 2, 2.3, 2.6, 3, 3.5, 4] #old
    #z = [ 0.7, 1 , 1.35, 1.7, 2, 2.3, 2.6, 3, 3.5, 4]# older
    print expt.parameters
    s=pylab.loadtxt('%s/recon_stats.txt'%outdir)
    xrecon=s[:-1,0]; yrecon=s[:-1,1]
    yrecon_d=s[:-1,2]; yrecon_u=s[:-1,3]
    yrecon_rms = s[:-1,4]
    yrecon_avr = s[:-1,5]
    yrecon_rms_down = yrecon_rms
    
    for i in range(len(yrecon_rms)):
       if yrecon_avr[i] < yrecon_rms[i]:
          yrecon_rms_down[i] = 10**(numpy.log10(yrecon_avr[i]) -0.01)
    
    #print xrecon
    #print yrecon
    #sys.exit()
    
    z =[0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.15, 2.35, 2.55, 2.85, 3.15, 3.5]#oldest ;)
    num = int(num)
    #num = 1.
    print num
    z_min,z_max, z_m = get_z(num)
    dl = get_dl(z_m)
    Vmax = get_Vmax(z_min,z_max)
    dsdl = get_dsdl(z_m,dl)

    z_m = (z_min + z_max)/2
    dl = get_dl(z_m)
    
    print z_min,z_max
    print expt.kind
    area = SURVEY_AREA*sqDeg2sr 
    area1 = SURVEY_AREA*sqDeg2sr
    print SURVEY_AREA
    
    SMIN  = get_sbins(numpy.power(10,Lmin),z_m,dl)*1e6
    SMAX  = get_sbins(numpy.power(10,Lmax),z_m,dl)*1e6
    sigma,fsigma = numpy.log10(get_Lbins([SURVEY_NOISE,SURVEY_NOISE*5],z_m,dl,'muJy'))#*(1.4/3)**(-.7))
    sigma2,fsigma2 = numpy.log10(get_Lbins([450,2400],z_m,dl,'muJy'))
    print z_m,dl, sigma,fsigma
    print SURVEY_NOISE,SURVEY_NOISE*5,SURVEY_AREA
    
    #area = 6672*sqDeg2sr #143165.15  49996.49 59876.8861135
    sbin3 = get_sbins(10**bins2,z_m,dl)*1e6
    sbin4 = get_sbins(10**bins3,z_m,dl)*1e6
    L_s   = numpy.log10(get_Lbins(ssbin,z_m,dl,'muJy'))
    #print expt.binsMedian
    #sys.exit()
    
    xreco = get_Lbins(xrecon,z_m,dl,'muJy')#*(1.4/3)**(-.7)
    yreco = yrecon
    
    
    
    print SMIN,SMAX
    #sys.exit()
    lin_yrecon = numpy.log10(yreco)
    print xrecon
    lin_xrecon = numpy.log10(xreco)
    
    
    bins2 = numpy.arange(15.,25.,0.4)
    
    bins5,rho_5,er5,rho_6,er6 = numpy.loadtxt('cos_data/skads_s1_3_LF.el',unpack=True)
 
      
    #novak 2017
    L_n   =[21.77,22.15,22.46,22.77,23.09,23.34]
    rho_n =[-2.85,-2.88,-3.12,-3.55,-4.05,-4.63]
    L_ner_u =[0.23,0.18,0.19, 0.2, 0.21, 0.28]
    L_ner_d =[1.1, 0.15, 0.14, 0.12, 0.12, 0.048]
    rho_ner=[0.09, 0.03,0.0355,0.059,0.1,0.234]
    
    #novak 2018
    L_n2   =[21.77,22.24,22.68,23.16,23.69,24.34, 24.74,25.56]
    rho_n2 =[-2.84,-2.90,-3.34,-4.00,-4.92,-4.92, -5.22,-5.07]	    
    L_ner_u2  =[0.23,0.27, 0.34, 0.37, 0.35, 0.21, 0.31, 0.03]
    L_ner_d2  =[1.0, 0.24, 0.17, 0.14, 0.16, 0.30, 0.20, 0.50]
    rho_ner_u2=[0.08,0.024,0.038,0.083,0.28,0.28 , 0.45, 0.34]
    rho_ner_d2=[0.07,0.023,0.035,0.070,0.25,0.25 , 0.37, 0.30]
    #l*= (1.4/3)**(-.7)
    #l2*= (1.4/3)**(-.7)
    #L*= (1.4/3)**(-.7) 
       
 
    #Lbins2 = numpy.log10(lf.get_Lbins(sbin2[sbin2>0],z_m,dl))
    print z_min,z_max
   
    bins  = numpy.arange(fsigma-0.,26.2,0.4)#bin +0.2
    #bins  = numpy.arange(fsigma -4.2,30,0.4) #bin
    #s,s_z = numpy.loadtxt('sdss/dr12s_s1_vol.txt', unpack=True, usecols=(-1,2))
   
    fig,ax = plt.subplots()
    
    plt.rc('lines', linewidth=2)
    
    if modelFamily in ['LFdpl','LFdpl_dpl']:
        faint =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
    elif modelFamily in ['LFpl','LFdpl_pl']:
        faint =lumfuncUtils.powerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lnorm_2)
    else:
        faint = lumfuncUtils.lognormpl(numpy.power(10,bins3), Lstar_2,Lslope_2,Lsigma,Lnorm_2)
        
    
  
    dm =0.#05
    hm=0.
    
    plt.errorbar(bins5,rho_6,yerr=er6,fmt='sr',label=r'$\rm{25 >\log_{10}(L) > %s}$'%(20.8),markersize=10)
    plt.errorbar(bins5,rho_5,yerr=er5,fmt='*b',label=r'$\rm{25 >\log_{10}(L) > %s}$'%(20.8),markersize=10)
    plt.plot(bins3+ hm ,numpy.log10(faint)-dm,'--b' ,label='faint-end')   
    plt.errorbar(lin_xrecon+ hm ,lin_yrecon-dm ,fmt='-k', label='MAP',linewidth=2)
    plt.errorbar(lin_xrecon+ hm,numpy.log10(yrecon_avr)-dm,fmt='-g', label = 'Average')
    #plt.errorbar(numpy.array(L_n),numpy.array(rho_n) -0.4,xerr=[L_ner_d,L_ner_u],yerr=rho_ner, fmt='^c', label='Novak et al 2017')
    #plt.errorbar(numpy.array(L_n2),numpy.array(rho_n2) -0.4,xerr=[L_ner_d2,L_ner_u2],yerr=[rho_ner_d2,rho_ner_u2], fmt='sk', label='Total RLF Novak et al 2018')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    
   
    plt.fill_between(lin_xrecon+hm,numpy.log10((yrecon_d)),numpy.log10(yrecon_u) -dm, color='k',alpha=0.2)
    
   
    plt.axvline(fsigma,color='k',linestyle='dashed')
    #plt.axvline(sigma,color='g')
    #plt.axvline(Llmin,color='r',linewidth=5)
   
    plt.xlabel(r'$\rm{log_{10}[L_{3 GHz}/(W}$ $\rm{Hz^{-1})]}$',fontsize = 30)
    plt.ylabel(r'$\rm{log_{10}[\rho_m(Mpc^{-3} mag^{-1})]}$',fontsize = 30)
    plt.text(19.9,-7.7,'%.1f < z < %.1f'%(0.1,0.4),fontsize = 30,alpha=0.6 ) # \n $M_i$ < -22
    #plt.text(23.9,-9.7,'%s'%outdir,fontsize = 16 ) # \n $M_i$ < -22
    plt.tick_params(axis='both',which = 'major', labelsize=20,width =3)
    plt.tick_params(axis='both',which = 'minor', labelsize=10, width=2)
    
    plt.xlim(18,26.)
    plt.ylim(-7,-1)
    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    #print truth['LMIN']
    plt.legend(frameon=False).draggable()
    plotf='%s/LF_recon_s1.pdf' % (outdir)
    plt.show()
    
	
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
