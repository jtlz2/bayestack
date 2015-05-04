#!/usr/bin/env python

"""
Simulate data set for bayestack.py
"""

import os,sys,shutil
import importlib
from utils import *
from countUtils import simulate

__name_cached=__name__
param_file=sys.argv[-1]

#setf='%s' % param_file.split('.')[-2]
# Import the settings variables
if __name__=='__main__':
    setf='bayestack_settings'
else:
    setf='%s'%param_file

set_module=importlib.import_module(setf)
globals().update(set_module.__dict__)
__name__=__name_cached

#try:
#    execfile(param_file)
#except IOError:
#    from settings import *

#-------------------------------------------------------------------------------

def main():
    """
    """

    global OUTPUT,DUMP

    if not os.path.exists(outdir): os.mkdir(outdir)
    OUTPUT=os.path.join(outdir,OUTPUT)
    if DUMP: DUMP=os.path.join(outdir,DUMP)

    shutil.copy(param_file,outdir)
    print 'Settings file: %s' % param_file

    r=simulate(simFamily,simParams,simParamsList,\
                          simBins,seed=SEED_SIM,N=NSIM,area=AREA_SIM,\
                          noise=NOISE_SIM,dump=DUMP,output=OUTPUT,\
                          verbose=True)
        
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

