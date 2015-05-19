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
    expt=countModel(modelFamily,nlaws,setf,dataset,binStyle,floatNoise)
    # Insert hacks here
    #plotRanges['C']=[0,200]
    labelDict=dict((name,name) for name in expt.parameters)
    plotTruth=dict((name,-99.0) for name in expt.parameters)
    plotRanges=dict((k,v[-2:]) for (k,v) in expt.priorsDict.items())

    # Load the data
    chain=pylab.loadtxt('%s/1-post_equal_weights.dat'%outdir)

    # Local plot
    line=True
    autoscale=False
    title='%s - %s'%(outdir,dataset)
    truth=plotTruth
    bundle=contour_plot.contourTri(chain,\
                            line=line,outfile='%s/%s'%(outdir,triangle),\
                            col=('red','blue'),labels=expt.parameters,\
                            ranges=plotRanges,truth=truth,\
                            autoscale=autoscale,title=title,\
                            labelDict=labelDict)

    # Plot for publication
    line=False
    autoscale=True
    title=''
    truth=None
    extn='pdf'
    bundle=contour_plot.contourTri(chain,\
                            line=line,\
                            outfile='%s/triangle-%s-publn.%s'%(outdir,run_num,extn),\
                            labels=expt.parameters,\
                            ranges=plotRanges,truth=truth,\
                            autoscale=autoscale,title=title,\
                            binsize=50,labelDict=labelDict)

    stats=fetchStats(outdir,expt.parameters,plotTruth)
    printLaTeX(expt.parameters,stats,dump=outdir)

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
