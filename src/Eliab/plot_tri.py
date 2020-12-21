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

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file
#param_dir=sys.argv[-1]
#param_file='settings.py'
#param_file=os.path.join(param_dir,param_file)
#try:
#    execfile(param_file)
#except IOError:
#    print 'IOError :('
#    #from settings import *

#print __name__
#__name_cached=__name__
#if __name__=='__main__':
#    param_file=sys.argv[-1]
#    settingsf=param_file.split('.')[-2]
#    set_module=importlib.import_module(settingsf)
#    globals().update(set_module.__dict__)
#    __name__=__name_cached

#-------------------------------------------------------------------------------


@profile
def main():

    """
    """

    print 'Settings file is %s' % setf

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    # Load prior ranges etc.
    expt=countModel(modelFamily,nlaws,setf,dataset,floatNoise)
    # Insert hacks here
    #plotRanges['C']=[0,200]
    print modelFamily
    if modelFamily=='LFsch':
    	labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$',\
    				'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha$'}#,'LZEVOL':'$Z_{evol}$'}
    
    elif modelFamily=='LFpl':
        print 'this..'
    	labelDict= {'LMIN':'$L_{min}$','LMAX':'$L_{max}$',\
    				'LNORM_2':'$\Phi_*$','LSTAR_2': '$L_*$','LSLOPE_2':r'$\alpha_1$','LSLOPE2_2':r'$\alpha_1$'}
    	
    elif modelFamily=='LFdpl':
        print 'this..'
    	labelDict= {'LMIN':'$L_{min}$','LMAX':'$L_{max}$',\
    				'LNORM_2':'$\Phi_*$','LSTAR_2': '$L_*$','LSLOPE_2':r'$\alpha_1$','LSLOPE2_2':r'$\alpha_2$','LZEVOL':'$Z_{evol}$'}
    elif modelFamily=='LFlognorm':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$',\
                                'LNORM_2':'$\Phi_*$','LSTAR_2': '$L_*$','LSLOPE_2':r'$\alpha_1$','LSIGMA':r'$\sigma$','LZEVOL':'$Z_{evol}$'}
                                
    elif modelFamily=='LFlognorm_dpl':
        print 'this is the right model'
        labelDict= {'LMIN':'$L_{min}$','LMAX2':'$L_{max2}$','LMIN2':'$L_{min2}$','LMAX':'$L_{max}$',\
                    'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\alpha_2$', 'LNORM_2':'$\Phi2_*$','LSTAR_2': '$L2_*$','LSLOPE_2':r'$\alpha_2$','LSIGMA':r'$\sigma$'}
        print labelDict
    elif modelFamily=='LFpl_lognorm':
        print 'this is the right model'
        labelDict= {'LMIN':'$L_{min}$','LMAX2':'$L_{max2}$','LMIN2':'$L_{min2}$','LMAX':'$L_{max}$',\
                    'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha_1$',\
                    'LNORM_2':'$\Phi2_*$','LSTAR_2': '$L2_*$','LSLOPE_2':r'$\alpha_2$','LSIGMA':r'$\sigma$'}
        print labelDict
    
    elif modelFamily=='LFdpl_pl':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX2':'$L_{max2}$','LMIN2':'$L_{min2}$','LMAX':'$L_{max}$',\
                    'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\alpha_2$', 'LNORM_2':'$\Phi2_*$','LSTAR_2': '$L2_*$','LSLOPE_2':r'$\beta_1$','LSLOPE2_2':r'$\beta_2$'}
    elif modelFamily=='LFpl_dpl':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX2':'$L_{max2}$','LMIN2':'$L_{min2}$','LMAX':'$L_{max}$',\
                    'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\alpha_2$', 'LNORM_2':'$\Phi2_*$','LSTAR_2': '$L2_*$','LSLOPE_2':r'$\beta_1$','LSLOPE2_2':r'$\beta_2$'}                    
    #elif modelFamily=='LFlognorm_dpl':
     #   labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$',\
      #              'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha_1$',\
       #             'LSLOPE2': r'$\alpha_2$', 'LNORM_2':'$\Phi2_*$','LSTAR_2': '$L2_*$','LSLOPE_2':r'$\alpha_2$','LSIGMA':r'$\sigma$'}
    elif modelFamily=='LFdpl_dpl':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX2':'$L_{max2}$','LMIN2':'$L_{min2}$','LMAX':'$L_{max}$',\
                    'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha_1$','LSLOPE2': r'$\alpha_2$', 'LNORM_2':'$\Phi2_*$',\
                    'LSTAR_2': '$L2_*$','LSLOPE_2':r'$\beta_1$','LSLOPE2_2':r'$\beta_2$'}   
    
    elif modelFamily=='LFdpl_dpl_z':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX2':'$L_{max2}$','LMIN2':'$L_{min2}$','LMAX':'$L_{max}$',\
                    'LNORM':'$\Phi_*$','LSTAR': '$L_*$','LSLOPE':r'$\alpha_1$','LSLOPE2': r'$\alpha_2$', 'LNORM_2':'$\Phi2_*$',\
                    'LSTAR_2': '$L2_*$','LSLOPE_2':r'$\beta_1$','LSLOPE2_2':r'$\beta_2$', 'A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$'}
    elif modelFamily=='LFevol':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$'}
    elif modelFamily=='LFevol_dpl':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$'}
    elif modelFamily=='LFevol_logn':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$'}
    elif modelFamily=='LFevol_phi_logn':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$','A_B':r'$\alpha^D_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$','B_D':r'$\beta^D_{SF}$', 'B_agn':r'$\beta_{AGN}$'} 
                    
    elif modelFamily=='LFevol_phi_logn_mat':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$','A_B':r'$\alpha^D_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'A_D':r'$\alpha^D_{SF}$','B_D':r'$\beta^D_{SF}$'}                
                    
    elif modelFamily=='LFevol_dpl_s':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$','B_SF':r'$\beta_{SF}$'}
        
    elif modelFamily=='LFevol_logn_s':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$','B_SF':r'$\beta_{SF}$'}
    elif modelFamily=='LFevol_logn_mat' or modelFamily=='LFevol_logn_el':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_agn':r'$L_{agn}$','A_SF':r'$\L_{SF}$'}
    
    elif modelFamily=='LFevol_logn_slope':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_agn':r'$L_{agn}$','A_SF':r'$\L_{SF}$','LSLOPE2_2':r'$\alpha_2$'}
    elif modelFamily=='LFevol_logn_sigma':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_agn':r'$L_{agn}$','A_SF':r'$\L_{SF}$','LSIGMA':r'$\sigma$'}
    elif modelFamily=='LFevol_logn_lmin':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_agn':r'$L_{agn}$','A_SF':r'$\L_{SF}$','LMIN2': '$L_{min2}$'}
    elif modelFamily=='LFevol_logn_lnorm':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_agn':r'$L_{agn}$','A_SF':r'$\L_{SF}$','LSTAR_2': '$L_*$','LNORM':'$\Phi_*$'} 
    elif modelFamily=='LFevol_logn_all':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$','LSTAR_2': '$L_*$','LNORM':'$\Phi_*$','LSIGMA':r'$\sigma$','LSLOPE2_2':r'$\alpha_2$'} 
    elif modelFamily=='LFevol_logn_all_L':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$','LSTAR_2': '$L_*$','LNORM':'$\Phi_*$','LSIGMA':r'$\sigma$','LSLOPE2_2':r'$\alpha_2$','LMIN_1':'$L_{min1}$', 'LMIN_2':'$L_{min2}$' ,'LMIN_3':'$L_{min3}$','LMIN_4':'$L_{min4}$' } 
    elif modelFamily=='LFevol_logn_L':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$', 'LMIN_1':'$L_{min1}$', 'LMIN_2':'$L_{min2}$' ,'LMIN_3':'$L_{min3}$' ,'LMIN_4':'$L_{min4}$' ,'LMIN_5':'$L_{min5}$' ,'LMIN_6':'$L_{min6}$' ,'LMIN_7':'$L_{min7}$' ,'LMIN_8':'$L_{min8}$' ,'LMIN_9':'$L_{min9}$'}        
        
    else:
    	labelDict=dict((name,lables) for name in expt.parameters)
    
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    if modelFamily in ['LFsch','LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFdpl_dpl_z','LFlognorm_dpl']:
        Lnorm=expt.parameters[expt.parameters.index('LNORM')]
        Lstar=expt.parameters[expt.parameters.index('LSTAR')]
        Lslope=expt.parameters[expt.parameters.index('LSLOPE')]
        #Lzevol=expt.parameters[expt.parameters.index('LZEVOL')]
    if modelFamily in ['LFdpl_dpl','LFlognorm_dpl']:
        Lslope2=expt.parameters[expt.parameters.index('LSLOPE2')]
    if modelFamily in ['LFlognorm','LFdpl','LFpl','LFdpl_pl','LFpl_dpl', 'LFlognorm_dpl','LFdpl_dpl','LFdpl_dpl_z']:
        Lnorm=expt.parameters[expt.parameters.index('LNORM_2')]
        Lstar=expt.parameters[expt.parameters.index('LSTAR_2')]
        Lslope=expt.parameters[expt.parameters.index('LSLOPE_2')]
    if modelFamily in ['LFdpl_dpl']:
        Lslope2_2=expt.parameters[expt.parameters.index('LSLOPE2_2')]
    
    if modelFamily in ['LFdpl_dpl_z']:
        alpha_agn=expt.parameters[expt.parameters.index('A_agn')]
        alpha_SF=expt.parameters[expt.parameters.index('A_SF')]
        beta_agn=expt.parameters[expt.parameters.index('B_agn')]
        beta_SF=expt.parameters[expt.parameters.index('B_SF')]

    
    noise = expt.parameters[expt.parameters.index('noise')]
    plotRanges=dict((k,v[-2:]) for (k,v) in expt.priorsDict.items())
    
    furniture={}
    
    # Load the data
    chain=pylab.loadtxt('%s/1-post_equal_weights.dat'%outdir) #
    
    
    print 'start'
    #print chain[0][1:]
    chain = numpy.array([c[1:] for c in chain])
    print chain[0]
    #sys.exit()
    #chain=pylab.loadtxt('%s/1-ev.dat'%outdir) 
    print expt.parameters,'oka'
    print modelFamily
    

    # Local plot
    line=False
    print labelDict
    autoscale=True
    print 'this is noise',noise
    title=''#r'$\rm{%s \mu Jy}$'%noise#%s - %s'%(outdir,dataset)
    stats=fetchStats(outdir,expt.parameters,plotTruth)
    printLaTeX(expt.parameters,stats,dump=outdir)
    truth=plotTruth
    bundle=contour_plot.contourTri(chain,\
                            line=line,\
                            labels=expt.parameters[1:],\
                            ranges=plotRanges,truth=truth,\
                            autoscale=autoscale,title=title,\
                            X_LABEL_OFFSET=-0.45,Y_LABEL_OFFSET=-0.35,FONTSIZE=15.,\
                            labelDict=labelDict\
                            ,furniture=furniture)

    # Plot for publication
    line=False
    autoscale=True
    title=''
    truth= {'noise': 150,'LMIN':21.7,'LMAX':24.6,'LNORM':numpy.log10(1e-7),'LSTAR': 23.1,'LSLOPE':1.,'LSLOPE2':1.}
    plotRanges={'noise': [148., 160],'LMIN': [0.2e22, 1.2e22],'LMAX':[1e23,1e25],'LNORM':[1e15,1e17],'LSTAR': [1e22,5e24],'LSLOPE':[0,6],'LSLOPE2':[0,6]}
    extn='png'
    bundle=contour_plot.contourTri(chain,\
                            line=line,\
                            outfile='%s/triangle-%s-publn.%s'%(outdir,run_num,extn),\
                            labels=expt.parameters,\
                            ranges=plotRanges,truth=truth,\
                            autoscale=autoscale,title=title,\
                            X_LABEL_OFFSET=-0.45,Y_LABEL_OFFSET=-0.35,FONTSIZE=15.,\
                            binsize=50,labelDict=labelDict,furniture=furniture)

    stats=fetchStats(outdir,expt.parameters,plotTruth)
    printLaTeX(expt.parameters,stats,dump=outdir)

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
