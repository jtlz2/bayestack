#!/usr/bin/env python

"""
This is reconstruct_lf_polt for the use of the plot_sfrd
It reconstruct the 95% region of the SFRD and stores them into
a file read by the plot_sfrd* 
June 20 
Eliab based on Jon's reconstruct.py

Usage:
Works for LFs (Sch,DPL,DDPL)
./reconstruct.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
from scipy import stats
import pylab as plt
import pymultinest
from bayestackClasses import countModel
from utils import sqDeg2sr,calculate_confidence,peak_confidence,fetchStats,calculate_confidence2
from lumfuncUtils import get_sbins,get_z,get_dl,get_sfr_q, get_sfrd_z, get_q,sfrd_Behroozi, sfrd_Madau

param_file=sys.argv[-1]
settingsf='%s.bayestack_settings' % param_file

#-------------------------------------------------------------------------------

def main():

    """
    """

    # Import the settings variables
    print 'Settings file is %s' % param_file

    # Import the settings variables
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)

    # Set up the experiment
    #bins=[12.679592, 31.84969509, 80.00281696, 200.9579904, 504.78364938, 1010, 3267.95919972]# 0.4 bins
    chains=['01a_3','01b_3','01c_3','01d_3','01e','01f','01g','01h','01i','01j'] #chains_1912 including peter
    chains=['01a_3','01b_2','01c_2','01d_2','01e_2','01f_2','01g','01h_2','01i_2','01j_3'] #chains_2001 DR4 remove phi2
    chains=['01a_8','01b_8','01c_8','01d_8','01e_8','01f_8','01g_7','01h_8','01i_8','01j_7'] #chains_2001 DR4 dpl_dpl
    chains=['01a_9','01b_9','01c_9','01d_9','01e_9','01f_9','01g_9','01h_9','01i_9','01j_9'] #chains_2001 DR4 dpl_dpl no mass selection
    chains=['02a_4_1','02a_4_2','02b_4','02c_4','02d_4','02e_4','02f_4','02g_4','02h_4','02i_4','02j_3']
    chains=['02a_z_6x_1','02a_z_6x_2','02b_z_7x','02c_z_7x','02d_z_7x','02e_z_7x','02f_z_7x','02g_z_7x','02h_z_7x','02i_z_7x','02j_z_7x' ]
    chains=['02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3','02z_31_3']
    chains=['02a_z_8x','02b_z_8x','02c_z_8x','02d_z_8x','02e_z_8x','02f_z_8x','02g_z_8x','02h_z_8x','02i_z_8x','02j_z_8x']#fixed model newz area
    chains=['02a','02b','02c','02d','02e','02f','02g','02h','02i','02j'] #lgnorm #chains_2001 DR4 newz new_area
    #chains=['02a_3','02b_3','02c_3','02d_3','02e_3','02f_3','02g_3','02h_3','02i_3','02j_3'] #mass <10
    #chains=['02a_z_6e','02b_z_6e','02c_z_6e','02d_z_6e','02e_z_6e','02e_z_6e','02g_z_6e','02h_z_6e','02i_z_6e','02j_z_6e']# New fixed 6e (using a)
    #chains=['02a_z_6f','02b_z_6f','02c_z_6f','02d_z_6f','02e_z_6f','02e_z_6f','02g_z_6f','02h_z_6f','02i_z_6f','02j_z_6f']# New fixed 6f (using a)
    #chains=['02a_z_6g','02b_z_6g','02c_z_6g','02d_z_6g','02e_z_6g','02e_z_6g','02g_z_6g','02h_z_6g','02i_z_6g','02j_z_6g']# New fixed 6f (using a)
    #chains=['01a_5','01b_5','01c_5','01d_5','01e_5','01f_5','01g_5','01h_5','01i_5','01j_5'] #stellar mass>10.2
    chains=['02a_5','02b_5','02c_5','02d_5','02e_5','02f_5','02g_5','02h_5','02i_5','02j'] #stellar mass>10
    chains=['02a_3','02b_3','02c_3','02d_3','02e_3','02f_3','02g_3','02h_3','02i_3','02j_3'] #stellar mass<10. logn
    #chains=['01a','01b','01c','01d','01e','01f','01g','01h','01i','01j'] #dpl #chains_2001 DR4 newz new_area
    #chains=['02a_z_5x','02b_z_5x','02c_z_5x','02d_z_5x','02e_z_5x','02f_z_5x','02g_z_5x','02h_z_5x','02i_z_5x','02j_z_5x']#stellar mass limit mass<10
    #chains=['02a_z_5f','02b_z_5f','02c_z_5f','02d_z_5f','02e_z_5f','02f_z_5f','02g_z_5f','02h_z_5f','02i_z_5f','02j_z_5x']#stellar mass limit mass<10
    #chains=['02a_z_8x_6_1','02b_z_8x_6_1','02c_z_8x_6_1','02d_z_8x_6_1','02e_z_8x_6_1','02f_z_8x_6_1','02g_z_8x_6_3','02h_z_8x_6_10','02i_z_8x_6_10','02j_z_8x_6_10']#fixed cut data more bins _el lmin lose
    #chains=['02a_z_8x_6_1','02b_z_8x_6_1','02c_z_8x_6_1','02d_z_8x_6_1','02e_z_8x_6_1','02f_z_8x_6_1','02g_z_8x_6_1','02h_z_8x_6_3','02i_z_8x_6_3','02j_z_8x_6_4']
    #chains='02z_31_5_15'
    pre_chain='11'
    bins=numpy.logspace(0,3.01,10)
    bins = list(bins)
    bins.append(2e3)
    bins1 = numpy.arange(15.2,30.2,0.1)
    nbins = len(bins)
    z_m = redshifts[0]
    dl = get_dl(z_m) 
    sbin1 = get_sbins(10**bins1,z_m,dl)*1e6
    expt=countModel(modelFamily,nlaws,settingsf,[dataset],floatNoise,doRedshiftSlices=True,mybins=sbin1)
    

    Lmin = 21#21.5# 19. 21.5
    Lmax=26.
    z_med=[]
    SFRD,SFRD_l,SFRD_h=[],[],[]
    sfrd_map=[]
    nov_z,nov_z_u, nov_z_l,  nov_sfrd, nov_sfrd_u, nov_sfrd_l,  Lower, lw_up, lw_lw = numpy.loadtxt('novak_2017.txt',usecols=(0,1,2, 3,4,5, 6,7,8),unpack=True)
    print modelFamily
    evaluations=[]
    evaluate=True
    for num in range(1,11):
     print num
     
     z_min,z_max, z_m = get_z(num,z_new=False)
     z_med.append(z_m)
     print z_m
    

     #f='%sev.dat'% outstem
     f='1-post_equal_weights.dat'
     f=os.path.join('chains_20%s%s'%(pre_chain,chains[num-1]),f)
     #f=os.path.join('chains_20%s%s'%(pre_chain,chains),f)
     print 'os.pathf',f
     
     #sfrd_map
     # Load equally-weighted posterior samples
     x=numpy.genfromtxt(f)
     nsamp=x.shape[0]
     ncols=x.shape[1] # The fifth [seventh] column is the posterior value
     # There must be a better way, but:
     #ncols = 14
     z=numpy.zeros((nsamp,ncols-1+expt.nbins))
     z[:,:-(expt.nbins-1)]=x
     # Shift posterior values to end
     z[:,-1]=z[:,ncols-1] # Copy...
     z[:,ncols-1]=0.0     # ...and blank
     
     print 'outstem',outstem

     # Fetch best-fit parameters and calculate best-fit line
     ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join('chains_20%s%s'%(pre_chain,chains[num-1]),outstem))
     drawmap=ana.get_best_fit()['parameters']
     
     #ana=pymultinest.analyse.Analyzer(ncols-1,\
     #   outputfiles_basename=os.path.join('chains_20%s%s'%(pre_chain,chains),outstem))
     #drawmap=ana.get_best_fit()['parameters']
    
     if True:
        print '--> Calculating *ML* reconstruction'
        drawmap=drawmap

     #ymap=expt.evaluate(drawmap)[0]#,expt.binsMedian)
     q_l,q,q_h = get_q(z_m)
     yrecon,lin_xrecon2,str_q_l=0,0,0
     sfrd=[]
     sfrd_map.append(get_sfrd_z(yrecon,lin_xrecon2,str_q_l,drawmap,expt.parameters,modelFamily,q_h ,26.,Lmin,z_m)[0])
     print drawmap
     for isamp in xrange(nsamp):
        #z[isamp,ncols-1:]=expt.evaluate(z[isamp,:])[0]#,expt.binsMedian)
        #Lmin =z[isamp,:][expt.parameters.index('LMIN')]

        Lmax =z[isamp,:][expt.parameters.index('LMAX')]
        #sfrd_z_l=get_sfrd_z(yrecon,lin_xrecon2,str_q_l,z[isamp,:],expt.parameters,modelFamily,q_l ,Lmax,Lmin,z_m)[0] #lower limit
        #sfrd_z_h=get_sfrd_z(yrecon,lin_xrecon2,str_q_l,z[isamp,:],expt.parameters,modelFamily,q_h ,Lmax,Lmin,z_m)[0] #upper limit
        sfrd_z  =get_sfrd_z(yrecon,lin_xrecon2,str_q_l,z[isamp,:],expt.parameters,modelFamily,q ,  Lmax,Lmin,z_m)[0]  #26, 21.5
        #plt.plot(z_m,numpy.log10(sfrd_z),'.', alpha=0.02)
        #plt.plot(z_m,numpy.log10(sfrd_z_l),'.', alpha=0.01)
        #plt.plot(z_m,numpy.log10(sfrd_z_h),'.', alpha=0.01)
        #sfrd.append(sfrd_z_l)
        sfrd.append(sfrd_z)
        #sfrd.append(sfrd_z_h)
     #plt.show()

    # Blanking, 0.0 -> NaN
    #sys.exit()

     s_med,dlow,dhigh,s_low,s_high=calculate_confidence(sfrd,alpha=0.95,ret_all=True)
     SFRD.append(s_med)
     SFRD_l.append(s_low)
     SFRD_h.append(s_high)
     '''
     plt.hist(numpy.log10(sfrd),20)
     plt.axvline(numpy.log10(s_med),color='r')
     plt.axvline(numpy.log10(s_low),color='k')
     plt.axvline(numpy.log10(s_high), color='k')
     plt.xlabel(r'$\rm{\log10(SFRD[M_\odot yr^{-1} Mpc^{-3}])}$',fontsize=20)
     plt.show()
     
    '''
    sfrdf='sfrd_chains_20%s%s_Lmin_%s'%(pre_chain,chains[0],'21_5')
    #sfrdf='sfrd_chains_20%s%s_Lmin_%s'%(pre_chain,chains,'21_5')
    fils=open(sfrdf,'w')
    hdr='# median redshift  sfrd_median  sfrd_low sfrd_high'
    fils.write('%s\n'%hdr)
    for i in range(len(SFRD)):
        fils.write('%f %f %f %f \n'%(z_med[i],SFRD[i],SFRD_l[i],SFRD_h[i]))
    fils.close()
    
    
    
    
    z_b=numpy.arange(0,60)/10.
    Ilbert_h = sfrd_Behroozi(z_b,z0=0.95-0.410,A=-1.,B=0.194-0.082,C=0.111+0.04)
    Ilbert_l = sfrd_Behroozi(z_b,z0=0.95+0.343,A=-1.,B=0.194+0.128,C=0.111-0.029)
    G_x , Gruppioni       = numpy.loadtxt('Gruppioni_2013',unpack=True,delimiter=',')
    G_x1, Gruppioni_yup   = numpy.loadtxt('Gruppioni_2013_y_up',unpack=True,delimiter=',')
    G_x2, Gruppioni_ydown = numpy.loadtxt('Gruppioni_2013_y_down',unpack=True,delimiter=',')
    Burg_z, Burg_FIR, Burg_FIR_er, Burg_tot, Burg_tot_er = numpy.loadtxt('Burgarella_2013',unpack=True)
    plt.fill_between(z_med,numpy.log10((SFRD_l)),numpy.log10(SFRD_h), color='b',alpha=0.2)
    plt.errorbar(z_med,numpy.log10(SFRD),fmt='*',markersize=14, label='Faint LF function')
    #plt.errorbar(z_med,numpy.log10(sfrd_map),fmt='pr',markersize=14, label='Faint LF function')
    plt.errorbar(z_b,numpy.log10(sfrd_Behroozi(z_b)),fmt='--',markersize=12, label='Behroozi et al. 2013')
    plt.errorbar(z_b,numpy.log10(sfrd_Behroozi(z_b,z0=0.95,A=-1.,B=0.194,C=0.111)),fmt=':',markersize=12, label='IIlbert et al. 2013 (UltraVista DR1)')
    plt.errorbar(z_b,numpy.log10(sfrd_Madau(z_b)),fmt='.-',markersize=12, label='Madau&Dickinson 2014')
    try:
        plt.errorbar(numpy.array(nov_z),numpy.array(nov_sfrd),xerr=[nov_z_l,nov_z_u],yerr=[nov_sfrd_l,nov_sfrd_u], fmt='sk', label='Total SFRD Novak et al 2017')
    except:
        print ''
    plt.ylabel(r'$\rm{\log10(SFRD[M_\odot yr^{-1} Mpc^{-3}])}$',fontsize=30)
    plt.xlabel(r'$\rm{z} $',fontsize = 30)    
    #plt.text(23.9,-9.7,'%s'%outdir,fontsize = 16 ) # \n $M_i$ < -22
    plt.tick_params(axis='both',which = 'major', labelsize=20,width =3)
    plt.tick_params(axis='both',which = 'minor', labelsize=10, width=2)
    plt.show()



    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
