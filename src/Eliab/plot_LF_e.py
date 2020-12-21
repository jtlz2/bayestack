#!/usr/bin/env python

"""
This is plot_counts.py
Jonathan Zwart
May 2014

Usage:

./plot_counts.py CHAINS_DIR

"""

import os,sys
#import os.path
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
from lumfuncUtils import get_Lbins,get_sbins,get_dl, get_Vmax,get_dsdl,get_z, get_dL
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

    expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise,doRedshiftSlices=True)
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
    chains=['a','a','b','c','d','e','f','g','h','i','j']
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
        #print datafiles, len(datafiles[0])
        dataf = datafiles[0]
        if len(dataf)==24:
            num = dataf[-5]
            print num
        else:
            num = dataf[-6:-4]
        d=numpy.genfromtxt('%s/%s'%(outdir,datafiles[0][8:]))
        #print '%s%s'%(outdir,datafiles[0][8:])   
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
    num=5
    print 'this is 1st num ', num
    #sys.exit()
    params = drawmap# [0]
    #for i in drawmap:
    #	params.append(i)
    print params
    print expt.parameters

    if modelFamily in ['LFsch','LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFdpl_dpl_z','LFlognorm_dpl','LFpl_lognorm']:
        Lnorm=params[expt.parameters.index('LNORM')]
        Lstar=params[expt.parameters.index('LSTAR')]
        Lslope=params[expt.parameters.index('LSLOPE')]
        #Lzevol=params[expt.parameters.index('LZEVOL')]
    if modelFamily in ['LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFdpl_dpl_z','LFlognorm_dpl']:
        Lslope2=params[expt.parameters.index('LSLOPE2')]
    if modelFamily in ['LFlognorm','LFpl', 'LFdpl', 'LFdpl_pl','LFpl_dpl', 'LFlognorm_dpl','LFdpl_dpl','LFdpl_dpl_z','LFpl_lognorm']:
        Lnorm_2=params[expt.parameters.index('LNORM_2')]
        Lstar_2=params[expt.parameters.index('LSTAR_2')]
        Lslope_2=params[expt.parameters.index('LSLOPE_2')]

    if modelFamily in ['LFdpl_dpl','LFdpl_dpl_z','LFdpl']:
        Lslope2_2=params[expt.parameters.index('LSLOPE2_2')]
    
    if modelFamily in ['LFlognorm','LFlognorm_dpl','LFpl_lognorm']:
        Lsigma = params[expt.parameters.index('LSIGMA')]
        
    if modelFamily in ['LFdpl_dpl_z','LFevol']:
	   alpha_agn=params[expt.parameters.index('A_agn')]
	   alpha_SF=params[expt.parameters.index('A_SF')]
	   beta_agn=params[expt.parameters.index('B_agn')]
	   beta_SF=params[expt.parameters.index('B_SF')]
    if modelFamily in ['LFevol_dpl_s','LFevol_logn_s']:
	   alpha_SF=params[expt.parameters.index('A_SF')]
	   beta_SF=params[expt.parameters.index('B_SF')]

    Lmin=params[expt.parameters.index('LMIN')]
    Lmax=params[expt.parameters.index('LMAX')]
    truth= {'noise': 150,'LMIN':22.,'LMAX':24.5,'LNORM':numpy.log10(1e-7),'LSTAR': 23.3,'LSLOPE':3.,'LSLOPE2':1.}    
    #z = [0, 0.5, 1 , 1.5, 2, 2.3, 2.6, 3, 3.5, 4] #old
    #z = [ 0.7, 1 , 1.35, 1.7, 2, 2.3, 2.6, 3, 3.5, 4]# older
    print expt.parameters
    s=pylab.loadtxt('%s/recon_stats_%s.txt'%(outdir,chains[num -1]))
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
    z_new=False
    Replot=False
    
    if z_new:
        z =[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,1.9,2.2,2.5, 2.8,3.2,3.6,4.0]
    else:
        z=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]#oldest ;)
    num = int(num)
    #num = 1.
    print num
    z_min,z_max, z_m = get_z(num,z_new=z_new)
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
    sigma,fsigma = numpy.log10(get_Lbins([SURVEY_NOISE,SURVEY_NOISE*5],z_m,dl,'muJy')*(1.4/3)**(-.7))
    print z_m,dl, sigma,fsigma
    print SURVEY_NOISE,SURVEY_NOISE*5,SURVEY_AREA
    
    #area = 6672*sqDeg2sr #143165.15  49996.49 59876.8861135
    sbin3 = get_sbins(10**bins2,z_m,dl)*1e6
    sbin4 = get_sbins(10**bins3,z_m,dl)*1e6
    #L_s   = numpy.log10(get_Lbins(ssbin,z_m,dl,'muJy'))
    #print expt.binsMedian
    #sys.exit()
    
    xreco = get_Lbins(xrecon,z_m,dl,'muJy')#*(1.4/3)**(-.7)
    yreco = yrecon
    
    
    
    print SMIN,SMAX
    #sys.exit()
    lin_yrecon = numpy.log10(yreco)
    lin_xrecon = numpy.log10(xreco)
    
    
    bins2 = numpy.arange(15.,25.,0.4)
    bins  = numpy.arange(fsigma,26.2,0.4)#bin +0.2
    if not os.path.isfile('cos_data/cos_s%s_LF.el'%num) and z_new or not os.path.isfile('cos_data/cos_s%s_LF_old.el'%num) and not z_new or Replot:
     print 'Reploting the RLF and storing them for next time.'
     #l,z_l = open_anytype('cos_data/cosmos_d1_0_8arc.txt',(12,21), F='uJy',L=True,getz=True)
     #l2,z_l2 = open_anytype('cos_data/cosmos_d1_0_8arc.txt',(12,24), F='uJy',L=True,getz=True)
     if z_new:
        L,z_L = open_anytype('cos_data/cos_s%s.txt'%num,(2,3), F='uJy',L=True,getz=True,band='S')
        print 'this is the new z bins (please check if lumfuncUtils and cos_manifest are using the right z)'
     else:
        L,z_L = open_anytype('cos_data/cos_s%s_old.txt'%num,(2,3), F='uJy',L=True,getz=True,band='S')
        print 'this is the old z binning style. (please check if lumfuncUtils and cos_manifest are using the right z)'
     L_c,z_c=open_anytype('cosmos_counter_parts_II.txt',(3,4),[z_min,z_max], F='uJy',L=True,getz=True, band='S')
     z_c2,S_c2,L_c2 =numpy.loadtxt('cosmos_counter_parts_II.txt',unpack=True, usecols=(3,4,5))
     L_c2=L_c2[(z_c2<z_max)&(z_c2 >z_min)]
     S_c2=S_c2[(z_c2<z_max)&(z_c2 >z_min)]
     z_c2=z_c2[(z_c2<z_max)&(z_c2 >z_min)]
    
    
     #rho_1,er1 = calc_Vmax(calcu_zeff(z_l,l,11.5e-32),l,bins,z_max,z_min,area)
     #rho_2,er2 = calc_Vmax(calcu_zeff(z_l2,l2,11.5e-32),l2,bins,z_max,z_min,area)
     #print numpy.log10(L)
     #print (1.4/3)**(-.7)
     #print numpy.log10(L*(1.4/3)**(-.7))
     #sys.exit()
     rho_3,er3 = calc_Vmax(calcu_zeff(z_L,L,1.15e-31),L,bins,z_max,z_min,area)
     rho_4,er4 = calc_Vmax(99.,L,bins,z_max,z_min,area)
     rho_5,er5 = calc_Vmax(calcu_zeff(z_c,L_c,1.15e-31),L_c,bins,z_max,z_min,area)
     rho_6,er6 = calc_Vmax(calcu_zeff(z_c2,10**L_c2,1.15e-31),10**L_c2,bins,z_max,z_min,area)
     if z_new:
        f = open('cos_data/cos_s%s_LF.el'%num,'w')
     else:
        f = open('cos_data/cos_s%s_LF_old.el'%num,'w')
     for i in range(len(bins)):
        f.write('%5.1f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n'%(bins[i], rho_3[i],er3[i],rho_5[i],er5[i],rho_6[i],er6[i]) )
     f.close()
    
     dl_c=get_dl(z_c2)
     L_c3=numpy.log10(get_Lbins(S_c2,z_c2,numpy.array(dl_c),'muJy')*(1.4/3.)**(-.7))
    else:
        if z_new:
            bins,rho_3,er3,rho_5,er5,rho_6,er6 = numpy.loadtxt('cos_data/cos_s%s_LF.el'%num, unpack=True)
            print 'Using stored RLF for the new z (The one with more z bins). If you want to calculate the plot the set "Replot=True"'
        else:
            bins,rho_3,er3,rho_5,er5,rho_6,er6 = numpy.loadtxt('cos_data/cos_s%s_LF_old.el'%num, unpack=True)  
            print 'Using stored RLF for the is old z bins. If you want to calculate the plot the set "Replot=True"'  
        
    
    
    
    
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
   
    
    #bins  = numpy.arange(fsigma -4.2,30,0.4) #bin
    #s,s_z = numpy.loadtxt('sdss/dr12s_s1_vol.txt', unpack=True, usecols=(-1,2))
    s,s_z = numpy.loadtxt('cos_data/cos_s%s.txt'%(num), unpack=True, usecols=(-1,2))
    #s=s*1e6
    print fsigma
    
    #sys.exit()
    
    Ncount = sum(Nbin[fsigma>ssbin])    
    Llmin = numpy.log10(get_Lbins([fsigma/numpy.sqrt(Ncount)],z_min,get_dl(z_min),'muJy'))[0]
    Llmin*=(1.4/3.)**(-.7)
    print Llmin
    
    L_s  = (get_Lbins(ssbin,z_m,dl,'muJy'))
    L_s*=(1.4/3.)**(-.7)
    L_s = numpy.log10(L_s)
    rho_1 = lumfuncs(Nbin,L_s ,calc_Vmax(calcu_zeff(z_m,10**L_s),L_s,L_s,z_max,z_min,area,getLF=False))
    #sys.exit()

   

    fig,ax = plt.subplots()
    
    plt.rc('lines', linewidth=2)
    '''
    if modelFamily in ['LFdpl','LFdpl_dpl']:
        faint =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
    elif modelFamily in ['LFpl','LFdpl_pl']:
        faint =lumfuncUtils.powerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lnorm_2)
    elif modelFamily in ['LFpl_dpl','LFdpl_dpl_z']:
        faint =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3),Lstar,Lslope,Lslope2,Lnorm)
    else:
        faint = lumfuncUtils.lognormpl(numpy.power(10,bins3), Lstar_2,Lslope_2,Lsigma,Lnorm_2)
    '''    
    
    #phi_dpl =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3),Lstar ,Lslope,Lslope2,Lnorm)
    #L_2 = numpy.power(10,bins3)/(1 + z_m)**(alpha_SF+z_m*beta_SF)
    L_2 = numpy.power(10,bins3)*(1.4/3.)**(-.7)/(1 + z_m)**(2.95-z_m*0.29)
    L_1 = numpy.power(10,bins3)*(1.4/3.)**(-.7)/(1 + z_m)**(2.86-z_m*0.7)
    phi_agn=lumfuncUtils.doublepowerlaw(L_1,24.29,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
    phi_sf =lumfuncUtils.doublepowerlaw(L_2,22.11 ,2.4,0.43,-3.10)#0.33 2.97
    phi_logn =lumfuncUtils.lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
    #phi_pl =lumfuncUtils.powerlaw(numpy.power(10,bins3),Lstar ,Lslope,Lnorm)
    dm =0.#05
    hm=0.
    Map =expt.evaluate(params)[0]
    bins_map = get_Lbins(expt.binsMedian,z_m,dl,'muJy')
    

    #plt.plot(numpy.log10(10**bins3 *(1.4/3)**(-.7)),numpy.log10(faint)-dm,'--b' ,label='faint-end')
    plt.plot(bins3,numpy.log10(phi_sf)-dm,'--c', label='local sf')
    plt.plot(bins3,numpy.log10(phi_agn)-dm,'--c', label='local agn')
    plt.plot(bins3,numpy.log10(phi_logn)-dm,'--r', label='local lognorm')  
    #plt.plot(numpy.log10(bins_map),numpy.log10(Map), label='my MAP')
    plt.errorbar(lin_xrecon+hm,lin_yrecon-dm ,fmt='-k', label='MAP',linewidth=2)
    plt.errorbar(lin_xrecon+hm,numpy.log10(yrecon_avr)-dm,fmt='-g', label = 'Average')
    #plt.plot(L_s2+.2, numpy.log10(rho_10), '*c',label='average LF',markersize=13)#bin + 0.2
    plt.errorbar(bins,rho_3 ,yerr=er3,fmt='ob',label='extracted',markersize=8) #bin
    #plt.errorbar(bins2[17:],rho_3[17:] ,yerr=er3[17:],fmt='ob',label='extracted',markersize=8) #bin
    #plt.errorbar(bins,rho_5 ,yerr=er5,fmt='*r',label='Total',markersize=8) #bin
    #plt.errorbar(bins,rho_6 ,yerr=er6,fmt='>g',label='Total LF',markersize=10) #bin
    plt.plot(L_s, numpy.log10(rho_1), '*m',label='low LF',markersize=10)
    #plt.errorbar(bins[1:],rho_2[1:] ,yerr=er2[1:],fmt='ok',label='dectected catalog',markersize=10) #bin
    #plt.errorbar(bins2[1:]-0.2,rho_3[1:],yerr=er3[1:],fmt='sr',label='noisy',markersize=10) #bin
    #plt.errorbar(bins2[1:]-0.2,rho_4[1:],yerr=er4[1:],fmt='vm',label='noisy lim',markersize=10) #bin
    #plt.errorbar(numpy.array(L_n),numpy.array(rho_n) -0.4,xerr=[L_ner_d,L_ner_u],yerr=rho_ner, fmt='^c', label='Novak et al 2017')
    #plt.errorbar(numpy.array(L_n2),numpy.array(rho_n2) -0.4,xerr=[L_ner_d2,L_ner_u2],yerr=[rho_ner_d2,rho_ner_u2], fmt='sk', label='Total RLF Novak et al 2018')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    
    #plt.plot(bins3+0.4,numpy.log10(logn)-0.4 ,'--k', label='faint-end')
    
    #plt.plot(bins3,numpy.log10(logn) ,'--', label='lognorm')
    
    #plt.plot(bins3,numpy.log10(logn2)-dm,'--r', label='lognorm2')
    #
    
    
    #plt.plot(bins3,numpy.log10(phi_dpl2), label='dpl bright-end2')
    #plt.plot(Lbins4[(Lbins4<23.6) &(Lbins4>22.12)],numpy.log10(phi_dpl)[(Lbins4<23.6) &(Lbins4>22.12)], '-c',linewidth=5 ,label = r'$\rm{ MAP dpl < 5\sigma}$')
    #plt.plot(bins2,numpy.log10(logn_2),label='starforming_2')
    #plt.plot(bins2,numpy.log10(logn_3),label='starforming_3')
    #plt.plot(bins2,numpy.log10(logn_4),label='starforming_4')
    
    #plt.fill_between(xrecon,numpy.log10(yrecon - yrecon_down2 ) ,numpy.log10(yrecon+yrecon_up2), color='c',alpha=0.2)
   
    #plt.errorbar(lin_xrecon+0.4,lin_yrecon -0.4,fmt='-b', label='MAP',linewidth=2)
    #plt.errorbar(lin_xrecon+0.4,lin_yrecon- 0.4,yerr=numpy.log10(yrecon_rms),fmt='r',label='rms')
    #plt.plot(L_s, numpy.log10(rho_1), '*c',label='average LF',markersize=10)
    #plt.plot(L_s5, numpy.log10(rho_14)-0.2, '*r',label='average LF',markersize=13)
    #plt.plot(L_s6+0.0, numpy.log10(rho_16)-0.2, '*b',label='average LF(bin)',markersize=13)
    #plt.plot(bins2[:-1], numpy.log10(rho_15), '*g',label='LF average',markersize=13)
    #print len(rho_11),len(rho_10),len(rho_12)
    #plt.fill_betweenx(numpy.log10(rho_10) ,L_s3,L_s4, color='red',alpha=0.2)
    
    
    #plt.plot(L_s6, numpy.log10(rho_13), '*b',label='average LF(bin+0.2)',markersize=13)#bin + 0.2
    #plt.plot(L_s1, numpy.log10(rho_1), 'oc',label='average LF(bin)',markersize=13)#bin
    
    #plt.plot(Lbins4,lin_phi ,'--b',label='MAP')
    
    #plt.fill_between(lin_xrecon+.0,numpy.log10((yrecon_avr -yrecon_rms))-dm,numpy.log10(yrecon_avr +yrecon_rms)-dm , color='r',alpha=0.3)
    
    #plt.fill_between(lin_xrecon[1:-1]+0.2,numpy.log10(yreco[1:-1] - yrecon_down[1:-1] )-dm,numpy.log10(yreco[1:-1]+yrecon_up[1:-1])-dm , color='b',alpha=0.2)
    plt.fill_between(lin_xrecon+hm,numpy.log10((yrecon_d)),numpy.log10(yrecon_u) -dm, color='k',alpha=0.2)
    
    #plt.errorbar(bins[2:],rho_1[2:],yerr=er1[2:],fmt = 'ob',label='Detected')
    #plt.errorbar(bins2,rho_6-0.2,yerr=er6,fmt='ob',label='Detected') #working one 
    #plt.plot(Lbins2,rho,'ob',label='simulation')
    #plt.plot(Lbins3,rho_,label='Noisy fluxes')
     #rho_2
    #plt.errorbar(bins,rho_2 -0.2,yerr=er2,fmt='ok',label='data',markersize=10) #bin +0.2
    
    
    #plt.errorbar(bins2,rho_6-0.2,yerr=er6,fmt='vm',label='Noisy NVSS') 
    #plt.errorbar(bins2,rho_7-0.2,yerr=er6,fmt='xc',label='NVSS Detected') 
    #for i in range(len(bins)):
    # print bins[i], rho_2[i]
    #sys.exit()
    
    #plt.errorbar(Lbins2,rho, err,fmt = 'ob',label='simulation')
    #plt.errorbar(Lbins2,rho_1,err_1,fmt = 'or',label='simulation_noisy')
    #plt.errorbar(L2,rho2,color='g',label='Analytic')
    #plt.errorbar(L3,rho3,color='g',label='Analytic II')
    #plt.axvspan(22.55 ,23.51,color='r',alpha=0.2)
    #plt.axvline(Lmin,color='r',linestyle='dashed')
    plt.axvline(fsigma,color='k',linestyle='dashed')
    #plt.axvline(sigma,color='g')
    plt.axvline(Llmin,color='r',linewidth=5)
    #plt.axvline(sigma2,color='m')
    #plt.axvline(fsigma2,color='m')
    #plt.axvline(23.89,color='b')
    #plt.axvline(23.71,color='c')
    #plt.axvline(23.52,color='r')
    #plt.axvline(truth['LMIN'],color='r', linestyle='dashed', label = 'True params')
    #plt.axvline(truth['LMAX'],color='r', linestyle='dashed')
    #plt.axvline(truth['LSTAR'],color='b')
    #plt.text(truth['LSTAR']-0.2,-9., '$L_*$',fontsize=18)    
    #plt.axvline(Lmax,color='r', label= 'MAP params')
    #plt.axvline(Lstar,color='b',linestyle='dashed')
    #plt.axvline(Lmax, color='r',linestyle='dashed')
    #plt.text(Lstar-0.12,-10.5, '$L_*$',fontsize=18)
    #plt.text(fsigma+0.1,-6.7, '$5\sigma$',fontsize=25)
    #plt.text(sigma+0.1,-6.7, '$\sigma$',fontsize=25)
    plt.xlabel(r'$\rm{log_{10}[L_{1.4 GHz}/(W}$ $\rm{Hz^{-1})]}$',fontsize = 30)
    plt.ylabel(r'$\rm{log_{10}[\rho_m(Mpc^{-3} mag^{-1})]}$',fontsize = 30)
    plt.text(19.9,-7.7,'%.1f < z < %.1f'%(z_min,z_max),fontsize = 30,alpha=0.6 ) # \n $M_i$ < -22
    #plt.text(23.9,-9.7,'%s'%outdir,fontsize = 16 ) # \n $M_i$ < -22
    plt.tick_params(axis='both',which = 'major', labelsize=20,width =3)
    plt.tick_params(axis='both',which = 'minor', labelsize=10, width=2)
    #plt.tight_layout()
    #plt.xtick(fontsize=15)
    #plt.ytick(fontsize=15)
    plt.xlim(18,27.)
    plt.ylim(-11.2,0.)
    #plt.ylim(-11,-5.5)
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
