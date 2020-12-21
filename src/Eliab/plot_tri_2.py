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
from pylab import*
from profile_support import profile
from utils import *
import contour_plot
from bayestackClasses import countModel
import corner

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file
corner_tri=True
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
    	labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_\rm{min}]$','LMAX':r'$\log_{10}[L_\rm{max}]$',\
    				'LNORM':'$\log_{10}[\Phi_1^*]$','LSTAR': '$\log_{10}[L_1^*]$','LSLOPE':r'$\alpha$'}#,'LZEVOL':'$Z_{evol}$'}
    
    elif modelFamily=='LFpl':
        print 'this..'
    	labelDict= {'LMIN':r'$L_{\rm{min}_2}$','LMAX':r'$L_{\rm{max}_2}$',\
    				'LNORM_2':'$\Phi_2^*$','LSTAR_2': '$L_2^*$','LSLOPE_2':r'$\alpha_1$','LSLOPE2_2':r'$\alpha_1$'}
    	
    elif modelFamily=='LFdpl':
        print 'this..'
    	labelDict= {'LMIN':r'$L_{\rm{min}_2}$','LMAX':r'$L_{\rm{max}_2}$',\
    				'LNORM_2':'$\Phi_2^*$','LSTAR_2': '$L_2^*$','LSLOPE_2':r'$\alpha_1$','LSLOPE2_2':r'$\beta_1$','LZEVOL':'$Z_{evol}$'}
    
    elif modelFamily=='LFlognorm':
        labelDict= {'LMIN':r'$L_{\rm{min}_2}$','LMAX':r'$L_{\rm{max}_2}$',\
                                'LNORM_2':'$\Phi_2^*$','LSTAR_2': '$L_2^*$','LSLOPE_2':r'$\alpha_1$','LSIGMA':r'$\sigma_{LF}$','LZEVOL':'$Z_{evol}$'}
                                
    elif modelFamily=='LFlognorm_dpl':
        print 'this is the right model'
        labelDict= {'LMIN':r'$\log_{10}[L_{\rm{min}_1}]$','LMAX2':r'$\log_{10}[L_{\rm{max}_1}]$','LMIN2':r'$\log_{10}[L_{\rm{min}_2}]$','LMAX':r'$\log_{10}[L_{\rm{max}_2}]$',\
                    'LNORM':'$\log_{10}[\Phi_1^*]$','LSTAR': '$\log_{10}[L_1^*]$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\beta_1$', 'LNORM_2':'$\log_{10}[\Phi_2^*]$','LSTAR_2': '$\log_{10}[L_2^*]$','LSLOPE_2':r'$\alpha_2$','LSIGMA':r'$\sigma_{LF}$'}
        print labelDict
    
    elif modelFamily=='LFdpl_pl':
        labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_{\rm{min}_1}]$','LMAX2': r'$\log_{10}[L_{\rm{max}_1}]$','LMIN2':r'$\log_{10}[L_{\rm{min}_2}]$','LMAX':r'$\log_{10}[L_{\rm{max}_2]$',\
                    'LNORM':'$\log_{10}[\Phi_1^*]$','LSTAR': '$\log_{10}[L_1^*]$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\beta_1$', 'LNORM_2':'$\log_{10}[\Phi_2^*]$','LSTAR_2': '$\log_{10}[L_2^*]$','LSLOPE_2':r'$\alpha_2$','LSLOPE2_2':r'$\beta_2$'}
                    
       
    elif modelFamily=='LFdpl_dpl':
        labelDict= {'noise': '$\sigma$','LMIN': r'$L_{\rm{min}_1}$','LMAX2': r'$L_{\rm{max}_1}$','LMIN2':r'$L_{\rm{min}_2}$','LMAX':r'$L_{\rm{max}_2}$',\
                    'LNORM':'$\Phi_1^*$','LSTAR': '$L_1^*$','LSLOPE':r'$\alpha_1$',\
                    'LSLOPE2': r'$\beta_1$', 'LNORM_2':'$\Phi_2^*$','LSTAR_2': '$L_2^*$','LSLOPE_2':r'$\alpha_2$','LSLOPE2_2':r'$\beta_2$'}

    elif modelFamily=='LFevol_logn_L':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$', 'A_agn':r'$\alpha_{AGN}$',\
                    'B_SF':r'$\beta_{SF}$', 'B_agn':r'$\beta_{AGN}$', 'LMIN_1':'$L_{min1}$', 'LMIN_2':'$L_{min2}$' ,'LMIN_3':'$L_{min3}$' ,'LMIN_4':'$L_{min4}$' ,'LMIN_5':'$L_{min5}$' ,'LMIN_6':'$L_{min6}$' ,'LMIN_7':'$L_{min7}$' ,'LMIN_8':'$L_{min8}$' ,'LMIN_9':'$L_{min9}$'} 

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
    elif modelFamily=='LFevol_logn_mat' or modelFamily=='LFevol_logn_el':
        labelDict= {'noise': '$\sigma$','LMIN':r'$\log_{10}[L_{\rm{min}}]$','LMAX':r'$\log_{10}[L_{\rm{max}}]$','A_agn':r'$\alpha_{AGN}$','A_SF':r'$\alpha_{SF}$'}        
                    
    elif modelFamily=='LFevol_dpl_s':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$','B_SF':r'$\beta_{SF}$'}
        
    elif modelFamily=='LFevol_logn_s':
        labelDict= {'noise': '$\sigma$','LMIN':'$S_{min} (\mu Jy)$','LMAX':'$L_{max}$','A_SF':r'$\alpha_{SF}$','B_SF':r'$\beta_{SF}$'}
    elif modelFamily=='LFevol_logn_mat':
        labelDict= {'noise': '$\sigma$','LMIN':'$L_{min}$','LMAX':'$L_{max}$','A_agn':r'$L_{agn}$','A_SF':r'$\L_{SF}$'}


    else:
    	labelDict=dict((name,lables) for name in expt.parameters)
    
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    if modelFamily in ['LFsch','LFdpl_pl','LFdpl_dpl','LFlognorm_dpl']:
        Lnorm=expt.parameters[expt.parameters.index('LNORM')]
        Lstar=expt.parameters[expt.parameters.index('LSTAR')]
        Lslope=expt.parameters[expt.parameters.index('LSLOPE')]
        #Lzevol=expt.parameters[expt.parameters.index('LZEVOL')]
    if modelFamily in ['LFdpl_dpl','LFlognorm_dpl']:
        Lslope2=expt.parameters[expt.parameters.index('LSLOPE2')]
    if modelFamily in ['LFlognorm','LFdpl','LFpl','LFdpl_pl', 'LFlognorm_dpl','LFdpl_dpl']:
        Lnorm=expt.parameters[expt.parameters.index('LNORM_2')]
        Lstar=expt.parameters[expt.parameters.index('LSTAR_2')]
        Lslope=expt.parameters[expt.parameters.index('LSLOPE_2')]
    if modelFamily in ['LFdpl_dpl']:
        Lslope2_2=expt.parameters[expt.parameters.index('LSLOPE2_2')]
    
    noise = expt.parameters[expt.parameters.index('noise')]
    plotRanges=dict((k,v[-2:]) for (k,v) in expt.priorsDict.items())
    '''
    furniture={'TRUNCATE_phi':False,'TRUNCATE_phi_limit':2.0e7,\
               'Lstar':Lnorm ,'FONTSIZE':4,\
               'ROTATION':60.0,'FIGSIZE':(8.27,11.69), 'DPI':400,\
               'AXIS_LABEL_OFFSET':-0.3,'LOG_BINS':Lnorm,\
               'PADDING':0.05}
    '''
    furniture={'TRUNCATE_Lstar':False,'TRUNCATE_Lstar_limit':2.0e7,\
               'FONTSIZE':4,\
               'ROTATION':60.0,'FIGSIZE':(8.27,11.69), 'DPI':400,\
               'AXIS_LABEL_OFFSET':-0.3}
    
    #furniture={}
    # Load the data
    chain=numpy.loadtxt('%s/1-post_equal_weights.dat'%outdir) #
    stats=fetchStats(outdir,expt.parameters,plotTruth)

    print 'start'
    #print chain[0][1:]
    chain = numpy.array([c[1:] for c in chain])
    #chain=numpy.delete(chain,-1,1)
    print expt.parameters
    #print expt.parameters.keys() 
    remove=['noise']#,'LMAX2', 'LMAX']
    #parameters= [ele for ele in expt.parameters if ele not in remove]
    #chain =numpy.delete(chain,[1],1)
    #print chain
    #print parameters
    #print chains
    #print chain2
    #sys.exit()
    '''
    print chain[0]
    #sys.exit()
    #chain=pylab.loadtxt('%s/1-ev.dat'%outdir) 
    print expt.parameters,'oka'
    print modelFamily
    
    labels=[labelDict[param] for param in expt.parameters[1:] ]
    rc('xtick',labelsize=24)
    rc('ytick',labelsize=24)
    
    if corner_tri:
        figure= corner.corner(chains,bins=50,labels=labels,color='#0057f6',label_kwargs={'fontsize':32},plot_contours=True,fill_contours=True,contour_kwargs={'linewidth':100})#contour_kwargs={'levels',0.6})
    show()
    sys.exit()
    
    
    plotRanges={'noise': [148., 160],'LMIN': [20, 26],'LMAX':[23,30] ,'LMIN2': [20, 28],'LMAX2':[23,30],\
                'LNORM':[-10,2],'LSTAR': [24,30],'LSLOPE':[-5,5],'LSLOPE2':[-5,5], 'LNORM_2':[-10,2],\
                'LSTAR_2': [20,26],'LSLOPE_2':[-2,2],'LSLOPE2_2':[-5,1]}
    # Local plot
    line=False
    print labelDict
    autoscale=True
    print 'this is noise',noise
    title=''#r'$\rm{%s \mu Jy}$'%noise#%s - %s'%(outdir,dataset)
    
    truth=plotTruth
    bundle=contour_plot.contourTri(chain,\
                            line=line,\
                            labels=expt.parameters[1:],\
                            ranges=plotRanges,truth=truth,\
                            autoscale=autoscale,title=title,\
                            X_LABEL_OFFSET=-0.8,Y_LABEL_OFFSET=-0.7 ,\
                            labelDict=labelDict\
                            ,furniture=furniture)
    '''
    # Plot for publication
    print 'latex stuff'
    printLaTeX(expt.parameters,stats,dump=outdir)
    line=False
    autoscale=True
    title=''
    truth= {'noise': 150,'LMIN':21.7,'LMAX':24.6,'LNORM_2':numpy.log10(1e-7),'LSTAR_2': 23.1,'LSLOPE_2':1.,'LSLOPE2_2':1.}
    
    truth={par:stats[par][0] for par in expt.parameters}
    print 'This is the truth\n \n'
    print truth
    
    
    #plotRanges={'noise': [148., 160],'LMIN': [20, 26],'LMAX':[26,30] ,'LMIN2': [20, 30],'LMAX2':[20,30],\
    #            'LNORM':[-10,-2],'LSTAR': [24,30],'LSLOPE':[-5,5],'LSLOPE2':[-5,5], 'LNORM_2':[-10,-2],\
    #            'LSTAR_2': [20,26],'LSLOPE_2':[-5,5],'LSLOPE2_2':[-5,5]}
    extn='png'
    plotRanges={}

    bundle=contour_plot.contourTri(chain,\
                            line=line,\
                            outfile='%s/triangle_%s.%s'%(outdir,outdir[-3:],extn),\
                            labels=expt.parameters[1:],\
                            ranges=plotRanges,truth=truth,\
                            autoscale=autoscale,title=title,\
                            X_LABEL_OFFSET=-0.55,Y_LABEL_OFFSET=-0.45,FONTSIZE=24.,\
                            binsize=50,labelDict=labelDict)
                            #X_LABEL_OFFSET=-0.45,Y_LABEL_OFFSET=-0.35,FONTSIZE=15.,\ dpl_dpl 
                            #X_LABEL_OFFSET=-0.45,Y_LABEL_OFFSET=-0.35,FONTSIZE=20.,\ model C

    stats=fetchStats(outdir,expt.parameters,plotTruth)
    print stats
    printLaTeX(expt.parameters,stats,dump=outdir)

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
