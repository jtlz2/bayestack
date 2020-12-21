#!/usr/bin/env python

"""
This is plot.py
Jonathan Zwart
April 2014

Usage:

./plot.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
if os.getenv('PBS_O_HOST') not in [None,'baltasar']:
    from matplotlib import pyplot as plt
import pylab
from profile_support import profile
from utils import *
import contour_plot
from bayestackClasses import countModel

#-------------------------------------------------------------------------------


@profile
def print_X(modelFamily,parameters,outdir):

    """
    """

    #print 'Settings file is %s' % setf

        #plotRanges['C']=[0,200]
    #print modelFamily
    if modelFamily=='LFsch':
    	labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_\rm{min}]$','LMAX':r'$\log_{10}[L_\rm{max}]$',\
    				'LNORM':'$\log_{10}[\Phi_1^*]$','LSTAR': '$\log_{10}[L_1^*]$','LSLOPE':r'$\alpha$'}#,'LZEVOL':'$Z_{evol}$'}
    
    elif modelFamily=='LFpl':
        print 'this..'
    	labelDict= {'LMIN':r'$\log_{10}[L_\rm{min}]$','LMAX':r'$\log_{10}[L_\rm{max}]$',\
    				'LNORM_2':'$\log_{10}[\Phi_1^*]$','LSTAR_2': '$\log_{10}[L_1^*]$','LSLOPE_2':r'$\alpha_1$','LSLOPE2_2':r'$\alpha_1$'}
    	
    elif modelFamily=='LFdpl':
        print 'this..'
    	labelDict= {'LMIN':r'$\log_{10}[L_\rm{min}]$','LMAX':r'$\log_{10}[L_\rm{max}]$',\
    				'LNORM_2':'$\log_{10}[\Phi_1^*]$','LSTAR_2': '$\log_{10}[L_1^*]$','LSLOPE_2':r'$\alpha_1$','LSLOPE2_2':r'$\beta_1$','LZEVOL':'$Z_{evol}$'}
    elif modelFamily=='LFlognorm':
        labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_\rm{min}]$','LMAX':r'$\log_{10}[L_\rm{max}]$',\
                                'LNORM_2':'$\log_{10}[\Phi_1^*]$','LSTAR_2': '$\log_{10}[L_1^*]$','LSLOPE_2':r'$\alpha_1$','LSIGMA':r'$\sigma_{LF}$','LZEVOL':'$Z_{evol}$'}
                                
    elif modelFamily=='LFlognorm_dpl':
        #print 'this is the right model'
        labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_{\rm{min_1}}]$','LMAX2':r'$\log_{10}[L_{\rm{max_1}}]$','LMIN2':r'$\log_{10}[L_{\rm{min_2}}]$','LMAX':r'$\log_{10}[L_{\rm{max_6}}]$',\
                    'LNORM':'$\log_{10}[\Phi_1^*]$','LSTAR': '$\log_{10}[L_1^*]$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\beta_1$', 'LNORM_2':'$\log_{10}[\Phi_2^*]$','LSTAR_2': '$\log_{10}[L_2^*]$','LSLOPE_2':r'$\alpha_2$','LSIGMA':r'$\sigma_{LF}$'}
        
    
    elif modelFamily=='LFdpl_pl':
        labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_{\rm{min_1}}]$','LMAX2':r'$\log_{10}[L_{\rm{max_6}}]$','LMIN2':r'$\log_{10}[L_{\rm{min2}}]$','LMAX':r'$\log_{10}[L_\rm{max2}]$',\
                    'LNORM':'$\log_{10}[\Phi_1^*]$','LSTAR': '$\log_{10}[L_1^*]$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\beta_1$', 'LNORM_2':'$\log_{10}[\Phi_2^*]$','LSTAR_2': '$\log_{10}[L_2^*]$','LSLOPE_2':r'$\alpha_2$','LSLOPE2_2':r'$\beta_2$'}
        
    elif modelFamily=='LFevol_logn_mat' or modelFamily=='LFevol_logn_el':
        labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_{\rm{min}}]$','LMAX':r'$\log_{10}[L_{\rm{max}}]$','A_agn':r'$\alpha_{AGN}$','A_SF':r'$\alpha_{SF}$'}        
       
    elif modelFamily=='LFdpl_dpl':
        labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_{\rm{min_1}}]$','LMAX2':r'$\log_{10}[L_{\rm{max_1}}]$','LMIN2':r'$\log_{10}[L_{\rm{min_2}}]$','LMAX':r'$\log_{10}[L_{\rm{max_6}}]$',\
                    'LNORM':'$\log_{10}[\Phi_1^*]$','LSTAR': '$\log_{10}[L_1^*]$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\beta_1$', 'LNORM_2':'$\log_{10}[\Phi_2^*]$','LSTAR_2': '$\log_{10}[L_2^*]$','LSLOPE_2':r'$\alpha_2$','LSLOPE2_2':r'$\beta_2$'}

    else:
    	labelDict=dict((name,lables) for name in parameters)
    
    plotTruth=dict((name,-99.0) for name in parameters)
 
   
    # Load the data
    chain=pylab.loadtxt('%s/1-post_equal_weights.dat'%outdir) #
    
    #print chain[0][1:]
    chain2=numpy.array([c[5:] for c in chain])
    chain = numpy.array([c[1:] for c in chain])
    
    # Local plot
    stats=fetchStats2(outdir,parameters)
    
    #printLaTeX2(parameters,stats)
    
    return stats,labelDict

#-------------------------------------------------------------------------------
'''
#chains=['101a','101b','101c','101d','101e','101f','101g','101h','101i','101j'] #dpl #chains_2001 DR4 newz new_area
#chains=['102a','102b','102c','102d','102e','102f','102g','102h','102i','102j'] #lgnorm #chains_2001 DR4 newz new_area
#chains=['202a_z_8x_6_1','202b_z_8x_6_1','202c_z_8x_6_1','202d_z_8x_6_1','202e_z_8x_6_1','202f_z_8x_6_1','202g_z_8x_6_1','202h_z_8x_6_1','202i_z_8x_6_1','202j_z_8x_6_1']

#mass conntributions
#chains=['101a_5','101b_5','101c_5','101d_5','101e_5','101f_5','101g_5','101h_5','101i_5','101j'] #Model A stellar mass>10
#chains=['102a_3','102b_3','102c_3','102d_3','102e_3','102f_3','102g_3','102h_3','102i_3','102i_3'] #stellar mass<10. logn
#chains=['202a_z_5x','202b_z_5x','202c_z_5x','202d_z_5x','202e_z_5x','202f_z_5x','202g_z_5x','202h_z_5x','202i_z_5x','202j_z_8x']#stellar mass limit mass<10
'''

chains=['101a','101b','101c','101d','101e','101f','101g','101h','101i','101j'] #dpl #chains_2001 DR4 newz new_area
#chains=['102a','102b','102c','102d','102e','102f','102g','102h','102i','102j'] #lgnorm #chains_2001 DR4 newz new_area
#chains=['202a_z_8x_6_1','202b_z_8x_6_1','202c_z_8x_6_1','202d_z_8x_6_1','202e_z_8x_6_1','202f_z_8x_6_1','202g_z_8x_6_1','202h_z_8x_6_1','202i_z_8x_6_1','202j_z_8x_6_1']

#mass conntributions
#chains=['102a_5','102b_5','102c_5','102d_5','102e_5','102f_5','102g_5','102h_5','102i_5','102j'] #Model A stellar mass>10
#chains=['102a_3','102b_3','102c_3','102d_3','102e_3','102f_3','102g_3','102h_3','102i_3','102i_3'] #stellar mass<10. logn
#chains=['202a_z_5x','202b_z_5x','202c_z_5x','202d_z_5x','202e_z_5x','202f_z_5x','202g_z_5x','202h_z_5x','202i_z_5x','202j_z_5x']#stellar mass limit mass<10
z=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]
statss=[]
params=[]
print
for num in range(1,11):
    #ana=pymultinest.analyse.Analyzer(ncols-1,\
    #    outputfiles_basename=os.path.join('chains_200%s'%chains[num-1],outstem))
    #drawmap=ana.get_best_fit()['parameters']
    
    #params = drawmap# [0]
    #print len(params)
    outdir='chains_201%s'%chains[num-1]
    f = open('%s/bayestack_settings.py'%outdir, 'r')
    mod = f.readlines()[31].split()[0]
    print mod
    if 'logn' in mod:
        if 'dpl' in mod:
            model='LFlognorm_dpl'
        elif '_el' in mod:
            model='LFevol_logn_el' 
        elif 'pl_' in mod:
            model='LFpl_lognorm'
        elif 'phi' in mod:
            model='LFevol_phi_logn_mat'
        elif 'mat' in mod:
            model='LFevol_logn_mat'
        else:
            model='LFlognorm'
    elif 'LFdpl' in mod:
        if '_dpl' in mod:
            model='LFdpl_dpl'
        elif '_pl' in mod:
            model='LFdpl_pl'
        else:
            model = 'LFdpl'
    
    else:
        model='LFpl'
    with open('chains_201%s/params.tex'%chains[num -1]) as f:
        parameters=[line.split(None,1)[0] for line in f]
    #print type(parameters)
    #print parameters
    
    #print parameters
    #sys.exit()
    stats,labelD=print_X(model,parameters,outdir)
    statss.append(stats)
    params.append(parameters)

print 'size of statss', len(statss)
print '%25s & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ & $%s < z < %s$ \\\\ \\hline \n'%('Parameter',z[0],z[1] ,z[1],z[2] ,z[2],z[3] ,z[3],z[4] ,z[4],z[5] ,z[5],z[6] ,z[6],z[7]  ,z[7],z[8],z[8],z[9],z[9],z[10]    ),
for i in range(1,13):
    p = parameters[i] 
    
    print '%25s &$%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ & $%6.2f_{-%3.2f}^{+%3.2f}$ \\\\ [3pt]'%\
          (labelD[p], statss[0][p][1], statss[0][p][1] - statss[0][p][2], statss[0][p][3] - statss[0][p][1] , \
                      statss[1][p][1], statss[1][p][1] - statss[1][p][2], statss[1][p][3] - statss[1][p][1] , \
                      statss[2][p][1], statss[2][p][1] - statss[2][p][2], statss[2][p][3] - statss[2][p][1],\
                      statss[3][p][1], statss[3][p][1] - statss[3][p][2], statss[3][p][3] - statss[3][p][1], \
                      statss[4][p][1], statss[4][p][1] - statss[4][p][2], statss[4][p][3] - statss[4][p][1],\
                      statss[5][p][1], statss[5][p][1] - statss[5][p][2], statss[5][p][3] - statss[5][p][1],\
                      statss[6][p][1], statss[6][p][1] - statss[6][p][2], statss[6][p][3] - statss[6][p][1],\
                      
                      statss[7][p][1], statss[7][p][1] - statss[7][p][2], statss[7][p][3] - statss[7][p][1],\
                      statss[8][p][1], statss[8][p][1] - statss[8][p][2], statss[8][p][3] - statss[8][p][1],\
                      statss[9][p][1], statss[9][p][1] - statss[9][p][2], statss[9][p][3] - statss[9][p][1])  

'''
vertical table
print r'$z_{\rm{min}}]$ & z_{\rm{max}}]$ &',r'$\log_{10}[L_{\rm{min}}]$ &',r'$\log_{10}[L_{\rm{max}}]$ &',r'$\log_{10}[L_{\rm{min}2}]$ &',r'$\log_{10}[L_{\rm{max}2}]$ &','$\log_{10}[\Phi_1^*]$ &','$\log_{10}[L_1^*]$ &',r'$\alpha_1$ &',r'$\beta_1$ &', '$\log_{10}[\Phi_2^*]$ &', '$\log_{10}[L_2^*]$ &',r'$\alpha_2$ &',r'$\beta_2$ \\'
for i in range(7):
    print '$%s $ &'%z[i], '%s$'%(z[i+1]),
    for param in params[i]:
        val=statss[i][param]
        #print '& $%6.1f_{%6.2f}^{%6.2f}$' % (round(val[0],1),round(val[2],2),round(val[3],2)),
        print '& $%6.1f_{-%3.2f}^{+%3.2f}$' % (round(val[0],1),round(val[0] - val[2],2),round(val[3] - val[0],2)),
    print '\\\\'
'''
