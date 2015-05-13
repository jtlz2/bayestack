#!/usr/bin/env python

"""
This is simulate.py
Jonathan Zwart
May 2015

Simulate data set for bayestack.py

Usage:

./simulate.py SETTINGS_FILE.py

"""

import os,sys,shutil
import importlib
from utils import *
from countUtils import simulate

param_file=sys.argv[-1]
setf='%s' % param_file.split('.')[-2]
print '%s is using %s' % (__name__,setf)

#-------------------------------------------------------------------------------

def main():
    """
    """

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    if not os.path.exists(outdir): os.mkdir(outdir)

    shutil.copy(param_file,outdir)
    print 'Settings file: %s' % param_file

    r=simulate(simFamily,simParams,simParamsList,\
                          simBins,seed=SEED_SIM,N=NSIM,area=AREA_SIM,\
                          noise=NOISE_SIM,dump=os.path.join(outdir,DUMP),\
                          output=os.path.join(outdir,OUTPUT),\
                          verbose=True,skadsf=skadsFile,\
                          simarrayf=simArrayFile,pole_posns=simPolePosns)
        
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

