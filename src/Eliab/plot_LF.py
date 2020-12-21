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
from lumfuncUtils import get_Lbins,get_sbins,get_dl, get_Vmax,get_dsdl,get_z, get_dL,dntoLF
from plat import open_anytype, calc_Vmax, calcu_zeff,lumfuncs
from matplotlib.ticker import AutoMinorLocator
from cycler import cycler
from scipy.interpolate import interp1d
z_c, c = numpy.loadtxt('completeness', unpack=1)
c_z= interp1d(z_c,c)

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
    if modelFamily in ['LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFlognorm_dpl']:
        Lslope2=params[expt.parameters.index('LSLOPE2')]
    if modelFamily in ['LFlognorm','LFpl', 'LFdpl', 'LFdpl_pl','LFpl_dpl', 'LFlognorm_dpl','LFdpl_dpl','LFdpl_dpl_z','LFpl_lognorm']:
        Lnorm_2=params[expt.parameters.index('LNORM_2')]
        Lstar_2=params[expt.parameters.index('LSTAR_2')]
        Lslope_2=params[expt.parameters.index('LSLOPE_2')]

    if modelFamily in ['LFdpl_dpl','LFdpl','LFdpl_dpl_z']:
        Lslope2_2=params[expt.parameters.index('LSLOPE2_2')]
    
    if modelFamily in ['LFlognorm','LFlognorm_dpl','LFpl_lognorm','LFevol_logn_sigma']:
        Lsigma = params[expt.parameters.index('LSIGMA')]
    if modelFamily in ['LFdpl_dpl_z','LFevol_dpl','LFevol_logn','LFevol_logn_all','LFevol_logn_all_L']:
	   alpha_agn=params[expt.parameters.index('A_agn')]
	   alpha_SF=params[expt.parameters.index('A_SF')]
	   beta_agn=params[expt.parameters.index('B_agn')]
	   beta_SF=params[expt.parameters.index('B_SF')]
    if modelFamily in ['LFevol_dpl_s','LFevol_logn_s']:
	   alpha_SF=params[expt.parameters.index('A_SF')]
	   beta_SF=params[expt.parameters.index('B_SF')]
    if modelFamily in ['LFevol_logn_lnorm','LFevol_logn_all','LFevol_logn_all_L']:
       	Lnorm=params[expt.parameters.index('LNORM')]
        Lstar=params[expt.parameters.index('LSTAR')]
        
	   
    if modelFamily in ['LFevol_logn_mat','LFevol_phi_logn_mat','LFevol_logn_el','LFevol_logn_sigma','LFevol_logn_lmin','LFevol_logn_lnorm']:
        alpha_SF=params[expt.parameters.index('A_SF')]
        alpha_agn=params[expt.parameters.index('A_agn')]
    if modelFamily in ['LFevol_logn_slope','LFevol_logn_all','LFevol_logn_all_L']:
        alpha_SF=params[expt.parameters.index('A_SF')]
        alpha_agn=params[expt.parameters.index('A_agn')]
        Lslope2_2=params[expt.parameters.index('LSLOPE2_2')]

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
    
    
     #novak 2018
    L_n2   =[[21.77,22.24,22.68,23.16,23.69,24.34, 24.74,25.56],\
             [22.30,22.61,22.96,23.38,23.80,24.10, 24.55,25.14],\
             [22.61,22.86,23.14,23.45,23.82,24.14, 24.40,24.71],\
             [22.85,23.16,23.69,24.24,24.80,25.31, 25.96,26.69],\
             [23.10,23.38,23.86,24.36,24.86,25.35, 25.94,26.36],\
             [23.32,23.57,23.94,24.32,24.67,25.06, 25.47,25.96],\
             [23.54,23.75,24.13,24.45,24.90,25.27, 25.74,26.10],\
             [23.73,23.99,24.57,25.10,25.68,26.18, 26.83,27.51],\
             [24.01,24.26,24.76,25.26,25.91,26.19, 27.45],\
             [24.30,24.56,24.80,25.13,25.37,25.80, 25.97,26.49]]
             
 
    L_ner_u2  =[[0.23,0.27, 0.34, 0.37, 0.35, 0.21, 0.31, 0.03],\
                [0.11,0.20, 0.26, 0.25, 0.24, 0.35, 0.31, 0.15],\
                [0.08,0.15, 0.20, 0.22, 0.17, 0.18, 0.24, 0.28],\
                [0.072,0.33,0.38, 0.40, 0.48, 0.42, 0.41, 0.28],\
                [0.080,0.30,0.32, 0.32, 0.32, 0.33, 0.24, 0.34],\
                [0.068,0.22,0.24, 0.25, 0.29, 0.30, 0.28, 0.21],\
                [0.067,0.22,0.21, 0.26, 0.18, 0.17,0.075, 0.11],\
                [0.093,0.36,0.31, 0.30, 0.25, 0.28, 0.16,0.028],\
                [0.076,0.34,0.36, 0.38, 0.24, 0.47, 0.27],\
                [0.097,0.14,0.20, 0.16, 0.23, 0.092,0.22,0.026]]
                
    
    L_ner_d2  =[[1.0, 0.24, 0.17, 0.14, 0.16, 0.30, 0.20, 0.50],\
               [0.29, 0.21, 0.15, 0.16, 0.17, 0.06, 0.099,0.29],\
               [0.24, 0.17, 0.12, 0.10, 0.15, 0.15, 0.08,0.074],\
               [0.16, 0.24, 0.19, 0.17, 0.15, 0.09, 0.16, 0.33],\
               [0.30, 0.19, 0.18, 0.18, 0.18, 0.17, 0.26, 0.18],\
               [0.14, 0.18, 0.15, 0.14, 0.10, 0.095,0.11, 0.21],\
               [0.19, 0.14, 0.15, 0.11, 0.19, 0.19, 0.29, 0.28],\
               [0.40, 0.17, 0.21, 0.23, 0.27, 0.25, 0.37, 0.53],\
               [0.20, 0.17, 0.16, 0.14, 0.27, 0.042, 0.27],\
               [0.22, 0.16, 0.10, 0.14,0.072, 0.21, 0.08, 0.30]]

    rho_n2 =[[-2.84,-2.90,-3.34,-4.00,-4.92,-4.92, -5.22,-5.07],\
             [-2.95,-3.17,-3.46,-4.24,-4.76,-5.41, -5.23,-5.44],\
             [-2.95,-3.11,-3.44,-3.81,-4.37,-4.48, -4.90,-5.14],\
             [-2.99,-3.24,-3.89,-4.54,-5.27,-5.33, -5.69,-5.89],\
             [-3.32,-3.51,-4.06,-4.74,-5.28,-5.43, -6.08,-5.70],\
             [-3.29,-3.54,-4.02,-4.45,-5.11,-5.56, -5.37,-6.07],\
             [-3.43,-3.61,-4.19,-4.56,-5.09,-5.34, -5.70,-5.75],\
             [-3.88,-3.98,-4.64,-5.33,-5.73,-6.27, -6.14,-6.77],\
             [-3.85,-4.28,-5.05,-5.60,-5.68,-6.44, -6.59],\
             [-4.43,-4.91,-5.46,-5.58,-5.89,-6.23, -6.26,-6.91]]
       
    rho_ner_u2=[[0.08,0.024,0.038,0.083,0.28,0.28 , 0.45, 0.34],\
                [0.048,0.027,0.036,0.089,0.17,0.45, 0.34, 0.45],\
                [0.062,0.024,0.032,0.048,0.096, 0.11,0.18,0.25],\
                [0.042,0.018,0.035,0.073,0.180, 0.20,0.34,0.45],\
                [0.045,0.019,0.033,0.072,0.15,  0.17,0.45,0.25],\
                [0.040,0.022,0.033,0.055,0.13,  0.22,0.17,0.45],\
                [0.061,0.020,0.035,0.055,0.099, 0.14,0.22,0.25],\
                [0.044,0.022,0.044,0.11, 0.16, 0.34, 0.28,0.76],\
                [0.048,0.025,0.058,0.13, 0.24, 0.34, 0.45],\
                [0.084,0.058,0.11 ,0.24, 0.20, 0.28, 0.28,0.76]]
                
    rho_ner_d2=[[0.07,0.023,0.035,0.070,0.25,0.25 , 0.37, 0.30],\
                [0.043,0.026,0.033,0.074,0.16,0.37, 0.30, 0.37],\
                [0.054,0.023,0.030,0.043,0.079,0.086,0.17,0.22],\
                [0.038,0.017,0.032,0.062,0.17,0.19, 0.30, 0.37],\
                [0.041,0.019,0.031,0.062,0.11,0.16, 0.37, 0.22],\
                [0.037,0.021,0.031,0.049,0.098,0.20,0.16, 0.37],\
                [0.054,0.019,0.032,0.049,0.081,0.11,0.20, 0.22],\
                [0.040,0.021,0.040,0.085,0.15, 0.30, 0.25,0.52],\
                [0.043,0.024,0.051,0.10, 0.15, 0.30, 0.37],\
                [0.070,0.051,0.087,0.16,0.19, 0.25, 0.25, 0.52]]
    #print xrecon
    #print yrecon
    #sys.exit()
    z_new=True
    Replot=False
    loglike=0.
    loglike_nov=0.
    lmin3=[ 19.2,  19.8,  20.5,  20.8,  21.6,  22.2,  22.5,  23. ,  22.9,  23. ]
    lmin3=[ 19.7,  20.3,  21. ,  21.3,  22.1,  22.7,  23. ,  23.5,  23.4,  23.5]
    lmin3=[ 20.2,  20.8,  21.5,  21.8,  22.6,  23.2,  23.5,  24. ,  23.9,  24. ]
    lmin3=[ 20.7,  21.3,  22. ,  22.3,  23.1,  23.7,  24. ,  24.5,  24.4,  24.5]
    lmin3=[ 21.7,  22.3,  23. ,  23.3,  24.1,  24.7,  25. ,  25.5,  25.4,  25.5]
    lmin3=[ 21.2,  21.8,  22.5,  22.8,  23.6,  24.2,  24.5,  25. ,  24.9,  25. ]
    
    sm=[0.05162624,  0.0651364 ,  0.17144366,  0.19570955,  0.73934002,    1.58958387]
    sm=[0.16325651,   0.20597938,   0.54215247,   0.61888793,         2.33799844,   5.02670556]
    sm=[0.51626243,   0.651364  ,   1.71443663,   1.95709547,  7.39340024,  15.89583871]
    sm=[1.63256513,    2.05979383,    5.42152466,    6.18887928, 23.37998442,   50.26705564]    
    sm=[16.32565134,    20.5979383 ,    54.21524659,    61.88879277,  233.79984425,   502.67055637]
    sm=[5.16,    6.51,   17.14,   19.57, 73.93,  158.95,  206.55]
    
    
    if not z_new:
        z =[0.1, 0.4, 0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]
    else:
        z=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]#oldest ;)
    num = int(num)
    #num = 1.
    print num
    z_min,z_max, z_m = get_z(num,z_new=False)
    dl = get_dl(z_m)
    Vmax = get_Vmax(z_min,z_max)
    dsdl = get_dsdl(z_m,dl)

    #z_m = (z_min + z_max)/2
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
    bins  = numpy.arange(fsigma -0.2,29.2,0.4)#bin +0.2
    #bins =[23.6,23.9,24.1,24.3,24.9,25.6,26.,26.7,28] #2< z<3.2
    #bins =[23.6, 23.9,24.1,24.4,24.6,24.9,25.2,25.6,26.,26.4,27.,27.5,28] #3.2< z<4
    bins=numpy.array(bins)
    if not os.path.isfile('cos_data/cos_s%s_LF.el'%num) and z_new or not os.path.isfile('cos_data/cos_s%s_LF_old.el'%num) and not z_new or Replot:
     print 'Reploting the RLF and storing them for next time.'
     #l,z_l = open_anytype('cos_data/cosmos_d1_0_8arc.txt',(12,21), F='uJy',L=True,getz=True)
     #l2,z_l2 = open_anytype('cos_data/cosmos_d1_0_8arc.txt',(12,24), F='uJy',L=True,getz=True)
     if z_new:
        L,z_L = open_anytype('cos_data/cos_s%s.txt'%num,(2,3), F='uJy',L=True,getz=True,band='S')
        print 'this is the new z bins (please check if lumfuncUtils and cos_manifest are using the right z)'
     else:
        L,z_L = open_anytype('cos_data/cos_s%s.txt'%num,(2,3), F='uJy',L=True,getz=True,band='S')
        print 'this is the old z binning style. (please check if lumfuncUtils and cos_manifest are using the right z)'
     L_c,z_c=open_anytype('cosmos_counter_parts_II.txt',(3,4),[z_min,z_max], F='uJy',L=True,getz=True, band='S')
     L_c2,z_c2=open_anytype('cosmos_comp_vla.el',(2,3),[z_min,z_max], F='uJy',L=True,getz=True, band='S')
     L_c3,z_c3=open_anytype('cosmos_vla.el',(2,3),[z_min,z_max], F='uJy',L=True,getz=True, band='S')
     L_c4,z_c4=open_anytype('/home/eliab/Documents/OC/Project/cosmos/Cosmos/dr4_cos2015_vla.tx',(-4,-1),[z_min,z_max], F='uJy',L=True,getz=True, band='S')
     L_c5,z_c5=open_anytype('/home/eliab/Documents/OC/Project/cosmos/Cosmos/dr4_cos2015_vla.tx',(3,-1),[z_min,z_max], F='uJy',L=True,getz=True, band='S')
     #z_c2,S_c2,L_c2 =numpy.loadtxt('cosmos_counter_parts_II.txt',unpack=True, usecols=(3,4,5))
     #L_c2=L_c2[(z_c2<z_max)&(z_c2 >z_min)]
     #S_c2=S_c2[(z_c2<z_max)&(z_c2 >z_min)]
     #z_c2=z_c2[(z_c2<z_max)&(z_c2 >z_min)]
    
    
     #rho_1,er1 = calc_Vmax(calcu_zeff(z_l,l,11.5e-32),l,bins,z_max,z_min,area)
     #rho_2,er2 = calc_Vmax(calcu_zeff(z_l2,l2,11.5e-32),l2,bins,z_max,z_min,area)
     #print numpy.log10(L)
     #print (1.4/3)**(-.7)
     #print numpy.log10(L*(1.4/3)**(-.7))
     #sys.exit()
     rho_3,er3,b3,n3 = calc_Vmax(calcu_zeff(z_L,L,1.15e-31),L,bins,z_max,z_min,area,table=True)
     #rho_4,er4 = calc_Vmax(99.,L,bins,z_max,z_min,area)
     print 'COS15 VLA'
     rho_5,er5,b5,n5 = calc_Vmax(calcu_zeff(z_c,L_c,1.15e-31),L_c,bins,z_max,z_min,1.77*sqDeg2sr,table=True)
     print 'mass lim DR4-vla'
     rho_6,er6,b6,n6 = calc_Vmax(calcu_zeff(z_c2,L_c2,1.15e-31),L_c2,bins,z_max,z_min,1.5*sqDeg2sr,table=True)
     print 'DR4-vla'
     rho_7,er7,b7,n7 = calc_Vmax(calcu_zeff(z_c3,L_c3,1.15e-31),L_c3,bins,z_max,z_min,1.5*sqDeg2sr,table=True)
     print 'cos-DR4-vla DR4 photo-z'
     rho_8,er8,b8,n8 = calc_Vmax(calcu_zeff(z_c4,L_c4,1.15e-31),L_c4,bins,z_max,z_min,1.5*sqDeg2sr,table=True)
     print 'cos-DR4-vla COS15 photo-z'
     rho_9,er9,b9,n9 = calc_Vmax(calcu_zeff(z_c5,L_c5,1.15e-31),L_c5,bins,z_max,z_min,1.5*sqDeg2sr,table=True)
     for i in range(len(b3)):
        print '%5.2f %7d %10d %13d %17d %11d'%(b3[i],n3[i],n7[i],n8[i],n9[i],n5[i])
     if z_new:
        f = open('cos_data/cos_s%s_LF_5.el'%num,'w')
        #_2 paper _5 stellar mass 10.2
        print 'saving data for future use in cos_data/cos_s%s_LF_5.el'%num
     else:
        f = open('cos_data/cos_s%s_LF_old.el'%num,'w')
        print 'saving data for future use in cos_data/cos_s%s_LF_old.el'%num
     for i in range(len(bins)):
        f.write('%5.1f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f %5.3f\n'%(bins[i], rho_3[i],er3[i],rho_5[i],er5[i],rho_6[i],er6[i],rho_7[i],er7[i],rho_8[i],er8[i],rho_9[i],er9[i]) )
     f.close()
    
    else:
        if z_new:
            bins,rho_3,er3,rho_5,er5,rho_6,er6 ,rho_7,er7 ,rho_8,er8 ,rho_9,er9 = numpy.loadtxt('cos_data/cos_s%s_LF_2.el'%num, unpack=True)
            print 'Using stored RLF for the new z (The one with more z bins ) from cos_data/cos_s%s_LF.el. If you want to calculate the plot the set "Replot=True"'%num
        else:
            bins,rho_3,er3,rho_5,er5,rho_6,er6 ,rho_7,er7 ,rho_8,er8 ,rho_9,er9 = numpy.loadtxt('cos_data/cos_s%s_LF_old.el'%num, unpack=True)  
            print 'Using stored RLF for the is old z bins in fr cos_data/cos_s%s_LF.el. If you want to calculate the plot the set "Replot=True"'%num
        
    
    
    
    
    #novak 2017
    '''
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
    '''
 
    #Lbins2 = numpy.log10(lf.get_Lbins(sbin2[sbin2>0],z_m,dl))
    print z_min,z_max
   
    
    #bins  = numpy.arange(fsigma -4.2,30,0.4) #bin
    #s,s_z = numpy.loadtxt('sdss/dr12s_s1_vol.txt', unpack=True, usecols=(-1,2))
    s,s_z = numpy.loadtxt('cos_data/cos_s%s.txt'%(num), unpack=True, usecols=(-1,2))
    #s=s*1e6
    print fsigma
    
    L_s  = (get_Lbins(ssbin,z_m,dl,'muJy'))*(1.4/3.)**(-.7)
    L_s = numpy.log10(L_s)
    rho_1 = lumfuncs(Nbin,L_s ,calc_Vmax(calcu_zeff(z_m,10**L_s),L_s,L_s,z_max,z_min,area,getLF=False))
    #L_s = numpy.log10(10**L_s*(1.4/3.)**(-.7))
    #sys.exit()
    
    rho_11 = dntoLF(L_s, Nbin, dsdl, Vmax)
    #dn     = expt.realise(params)
    #nov_real=expt.realise([2.3,2.95,2.86,-0.29,-0.7,19.5,26.0])
    

    fig,ax = plt.subplots()
    
    plt.rc('lines', linewidth=2)
    
    if modelFamily in ['LFdpl','LFdpl_dpl','LFdpl_dpl_z']:
        faint =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
    elif modelFamily in ['LFpl','LFdpl_pl']:
        faint =lumfuncUtils.powerlaw(numpy.power(10,bins3), Lstar_2,Lslope_2,Lnorm_2)
    elif modelFamily in ['LFpl_dpl','LFdpl_dpl_z']:
        faint =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3),Lstar,Lslope,Lslope2,Lnorm)
    elif modelFamily =='LFevol_logn_mat':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_mat',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for lL in bins3]
    elif modelFamily =='LFevol_phi_logn_mat':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_phi_logn_mat',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for lL in bins3]
    elif modelFamily =='LFevol_logn_el':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_el',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for lL in bins3]
    elif modelFamily =='LFevol_phi_logn':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_phi_logn',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]
    elif modelFamily =='LFevol_logn':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]

    elif modelFamily =='LFevol_logn_slope':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_slope',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]
    elif modelFamily =='LFevol_logn_sigma':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_sigma',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]
    elif modelFamily =='LFevol_logn_lstar':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_lstar',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]
    elif modelFamily =='LFevol_logn_lmin':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_lmin',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]
    elif modelFamily =='LFevol_logn_lnorm':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_lnorm',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]
    elif modelFamily =='LFevol_logn_all':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_all',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]

    elif modelFamily =='LFevol_logn_all_L':
        faint =[lumfuncUtils.LF(L=numpy.power(10,lL),params=params,paramsList=expt.parameters,\
               family='LFevol_logn_all_L',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True,SF=True) for  lL in bins3]

    else:
        faint = lumfuncUtils.lognormpl(numpy.power(10,bins3), Lstar_2,Lslope_2,Lsigma,Lnorm_2)
    
    
    #phi_dpl =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3),Lstar ,Lslope,Lslope2,Lnorm)
    
    L_1 = numpy.power(10,bins3)/(1 + z_m)**(2.86-z_m*0.7)
    L_2 = numpy.power(10,bins3)/(1 + z_m)**(2.95-z_m*0.29)
    #L_4 = numpy.power(10,bins3)/(1 + z_m)**(alpha_SF)
    #L_3 = numpy.power(10,bins3)/(1 + z_m)**(alpha_agn)
    #L_4 = numpy.power(10,bins3)/(10)**(alpha_SF)
    #L_3 = numpy.power(10,bins3)/(10)**(1)#alpha_agn)
    #for i in range(len(L_3)):
        #print bins3[i],numpy.log10(L_3[i])
    #L_4 = numpy.power(10,bins3)/(1 + z_m)**(alpha_SF+z_m*beta_SF)
    #L_3 = numpy.power(10,bins3)/(1 + z_m)**(alpha_agn+z_m*beta_agn)
    phi_agn_n=lumfuncUtils.doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(10**(-5.5)))
    phi_sf_n =lumfuncUtils.lognormpl(L_2, numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3/2.5))#0.33 2.97   
    
    #phi_sf =lumfuncUtils.lognormpl(L_4, 20.56,1.22,0.63,-2.26)#numpy.log10(3.55e-3/2.5))#0.33 2.97 
    phi_sf2 =lumfuncUtils.lognormpl(L_2, numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3/2.5))#0.33 2.97 
    
    #phi_sf =lumfuncUtils.lognormpl(L_4, Lstar,1.22,0.63,Lnorm)#numpy.log10(3.55e-3/2.5))#0.33 2.97
    #phi_sf =lumfuncUtils.lognormpl(L_4, numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3/2.5))#0.33 2.97 
     

    
    nov_Tot_n= phi_agn_n+phi_sf_n 
    #phi_agn=lumfuncUtils.doublepowerlaw(L_3,24.59,1.27,0.49,numpy.log10(10**(-5.5)))
    #print 'this is Lminz', numpy.log10(10**18/(1 + z_m)**(alpha_SF+z_m*beta_SF))
    #phi_agn=lumfuncUtils.doublepowerlaw(L_3,24.59,1.27,0.49,numpy.log10(10**(-5.5)))
    phi_1 =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3),23.68 ,3.68 ,0.8, -3.36)#0.33 2.97    
    phi_2 =lumfuncUtils.doublepowerlaw(numpy.power(10,bins3), 24.5,1.37,0.5 ,-4.)
    #tot = phi_1+phi_2
    dm =0.4#02
    hm=0.
    #Map =expt.evaluate(params)[0]
    bins_map = get_Lbins(expt.binsMedian,z_m,dl,'muJy')
    bins1_4 = numpy.log10(10**bins*(1.4/3.)**(-.7))
    #fix=[lumfuncUtils.LF(L=numpy.power(10,lL),params=[3.75,2.55,4.74,0,0,20.3,26.65],paramsList=expt.parameters,\
    #           family='LFevol_logn',inta=None,area=0,S=0, z=z_m, dlds=0,Vmax=0,dl=0,bypassLim=True) for  lL in bins3]
    

    #plt.plot(numpy.log10(10**bins3 *(1.4/3)**(-.7)),numpy.log10(faint)-dm,'--b' ,label='faint-end')
    
    #plt.plot(numpy.log10(bins_map),numpy.log10(Map), label='my MAP')
    #plt.plot(bins3,numpy.log10(fix)-dm,'--c', label='Fix')
    #plt.plot(bins3,numpy.log10(phi_agn),'-.c', label='local agn')
    
    
    #plt.plot(bins3,numpy.log10(phi_sf),'--r', label='local sf')
    #plt.plot(bins3,numpy.log10(phi_sf2),'--c', label='local sf2')
    #plt.plot(numpy.log10(10**bins3),numpy.log10(phi_dpl)-dm,'--c', label='agn')
    #plt.plot(numpy.log10(10**bins3*(1.4/3.)**(.7)),numpy.log10(tot)-0.4,'-.g',linewidth=5, label='Tot')
    
    
    #plt.plot(bins3,numpy.log10(phi_agn),'--r', label='local agn nov')
    plt.errorbar(lin_xrecon+hm,lin_yrecon-dm,fmt='-b', label='This work Total',linewidth=2)
    plt.plot(numpy.log10(10**bins3),numpy.log10(faint)-dm,'--b', label='This work SFG')
    #plt.errorbar(lin_xrecon+hm,numpy.log10(yrecon_avr)-dm,fmt='-g', label = 'Average')
    #plt.plot(L_s2+.2, numpy.log10(rho_10), '*c',label='average LF',markersize=13)#bin + 0.2
    plt.errorbar(bins,rho_3 ,yerr=er3,fmt='ob',label='This work 1/Vmax',markersize=8) #bin
    #plt.errorbar(bins2[17:],rho_3[17:] ,yerr=er3[17:],fmt='ob',label='extracted',markersize=8) #bin
    #plt.errorbar(bins,rho_5 ,yerr=er5,fmt='^k',label='VLA-COSMOS2015',markersize=10) #bin
    #plt.errorbar(bins,rho_8 ,yerr=er8,fmt='sc',label='VLA-COSMOS2015 in DR4 using DR4 z',markersize=10) #bin
    #plt.errorbar(bins,rho_9 ,yerr=er9,fmt='*m',label='VLA-COSMOS2015 in DR4 using COS15z',markersize=14) #bin
    #plt.errorbar(bins,rho_6 ,yerr=er6,fmt='ok',label='DR4 detected',fillstyle='none',markersize=14) #bin
    #plt.errorbar(bins,rho_7 ,yerr=er6,fmt='dg',label='DR4 detected - Full NIR sample',markersize=12) #bin
    #plt.plot(L_s, numpy.log10(rho_1), '*m',label='low LF',markersize=10)
    #plt.plot(L_s, numpy.log10(rho_11), 'sc',label='reconstruct LF',markersize=10)
    #plt.errorbar(numpy.array(L_n2[num-1]),numpy.array(rho_n2[num-1]) -0.4,xerr=[L_ner_d2[num-1],L_ner_u2[num-1]],yerr=[rho_ner_d2[num-1],rho_ner_u2[num-1]], fmt='sr',fillstyle='none',ecolor='r',markeredgewidth=1.5, label='Total RLF Novak+2018')
    plt.errorbar(numpy.log10(10**bins3), numpy.log10(nov_Tot_n),fmt='-r',label='Total Novak + 2018')
    plt.plot(bins3,numpy.log10(phi_sf_n),'-.r', label='SFG Novak + 2017')
    '''
    print 'This is the dn'
    
    print dn
    
    for key in dn.keys():
        if round(key,2)==round(z_m,2):z_m=key

    dn=dn[z_m]
    nov_real=nov_real[z_m]
    for i in range(len(L_s)):
            loglike_i= Nbin[i]*numpy.log(dn[i]) + Nbin[i] - (Nbin[i] + 0.5)*numpy.log(Nbin[i]) - 0.5*numpy.log(2*numpy.pi)  - dn[i]
            loglike_nov_i= Nbin[i]*numpy.log(nov_real[i]) + Nbin[i] - (Nbin[i] + 0.5)*numpy.log(Nbin[i]) - 0.5*numpy.log(2*numpy.pi)  - nov_real[i]
            if nov_real[i]>0:loglike_nov+=loglike_nov_i
            if dn[i]>0:loglike+=loglike_i
            print '%10.5f %8.2d %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f'%(ssbin[i], Nbin[i], dn[i],nov_real[i],loglike_i,loglike_nov_i,loglike, loglike_nov)
    '''        
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
    plt.fill_between(lin_xrecon+hm,numpy.log10((yrecon_d))-dm,numpy.log10(yrecon_u) -dm, color='k',alpha=0.2)
    
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
    #plt.axvline(Llmin,color='r',linewidth=5)
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
    #plt.axvline(Lmin,color='b')
    #plt.axvline(Lmax, color='r',linestyle='dashed')
    #plt.text(21,-1, '$L_*$',fontsize=18)
    #plt.text(fsigma+0.1,-6.7, '$5\sigma$',fontsize=25)
    #plt.text(sigma+0.1,-6.7, '$\sigma$',fontsize=25)
    plt.xlabel(r'$\rm{log_{10}[L_{1.4 GHz}/(W}$ $\rm{Hz^{-1})]}$',fontsize = 30)
    plt.ylabel(r'$\rm{log_{10}[\rho_m(Mpc^{-3} mag^{-1})]}$',fontsize = 30)
    lmin3_14 = 10**lmin3[num-1]*(1.4/3.)**(-.7)
    #plt.text(22,-7.7,'%.1f < z < %.1f %s Lmin=%6.2f %s Lmin_cut=%6.2f %s Smin_cut=%6.2f $\mu$Jy'%(z_min,z_max,'\n',Lmin,'\n',numpy.log10(lmin3_14),'\n',sm[num-1]),fontsize = 30,alpha=0.6 ) # \n $M_i$ < -22
    #plt.text(22,-7.7,'%.1f < z < %.1f %s Lmin=%6.2f'%(z_min,z_max,'\n',Lmin),fontsize = 30,alpha=0.6 ) # \n $M_i$ < -22
    plt.text(22,-7.7,'%.1f < z < %.1f'%(z_min,z_max),fontsize = 30,alpha=0.6 ) # \n $M_i$ < -22
    #plt.text(22.9,-8.7,'%s'%outdir,fontsize = 16 ) # \n $M_i$ < -22
    plt.tick_params(axis='both',which = 'major', labelsize=20,width =3)
    plt.tick_params(axis='both',which = 'minor', labelsize=10, width=2)
    #plt.tight_layout()
    #plt.xtick(fontsize=15)
    #plt.ytick(fontsize=15)
    plt.xlim(20,28.)
    plt.ylim(-9,-2)
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
