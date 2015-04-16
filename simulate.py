#!/usr/bin/env python
##!/Library/Frameworks/EPD64.framework/Versions/Current/bin/python
##!/usr/bin/python

"""
Simulate data set for lumfunc.py
Output to XXXX.txt
"""

import os,sys,shutil
import importlib
from utils import *
import lumfunc


__name_cached=__name__
param_file=sys.argv[-1]

#setf='%s' % param_file.split('.')[-2]
# Import the settings variables
if __name__=='__main__':
    setf='settings'
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

    #if not simulate:
    #    print 'Simulation not requested!'
    #    return 1

    if not os.path.exists(outdir): os.mkdir(outdir)
    OUTPUT=os.path.join(outdir,OUTPUT)
    if DUMP: DUMP=os.path.join(outdir,DUMP)

    shutil.copy(param_file,outdir)
    note='Settings file: %s' % param_file
    print note

    if not batch_sim:
        dn_by_ds=lumfunc.simtable(bins,a=ALPHA_SIM,N=N_SIM,A=AREA_SIM,\
                        Smin=SMIN_SIM,Smax=SMAX_SIM,\
                        seed=SEED_SIM,noise=NOISE_SIM,dump=DUMP,output=OUTPUT,\
                        verbose=verbose,version=2,NLAWS_SIM=NLAWS_SIM,\
                        b=BETA_SIM,S0=S0_SIM)
        print 'Look in %s' % OUTPUT

    else:
        median_bins=medianArray(bins)
        simf='sims.dat'
        simf=os.path.join(outdir,simf)
        s=open(simf,'w')
        head1='# %s\n' % (' '.join(['%i'%(i+1) for i in range(nbins-1)]))
        s.write(head1)
        #head2='%s\n' % (' '.join(['%f'%i for i in median_bins]))
        #s.write(head2)
        refresh=50
        for ibatch in range(nbatch):
            iseed=SEED_SIM+ibatch
            idump='%s_%i.%s' % (DUMP.split('.')[0],iseed,DUMP.split('.')[-1])
            ioutput='%s_%i.%s' % (OUTPUT.split('.')[0],iseed,OUTPUT.split('.')[-1])
            idump=False; ioutput=None
            dn=lumfunc.simtable(bins,a=ALPHA_SIM,N=N_SIM,A=AREA_SIM,\
                        Smin=SMIN_SIM,Smax=SMAX_SIM,\
                        seed=iseed,noise=NOISE_SIM,dump=idump,output=ioutput,\
                        verbose=verbose,version=2)
            line='%s\n' % (' '.join(['%i'%i for i in dn]))
            s.write(line)
            if ibatch % refresh == 0: s.flush()
            print 'sim %i written to %s' % (ibatch,simf)

        s.close()

        
    return 0


#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

