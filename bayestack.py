#!/usr/bin/env python


"""
"""

import os
from bayestack_settings import *
from bayestackClasses import *
import pymultinest

#-------------------------------------------------------------------------------

def myprior(cube,ndim,nparams):
    cube=model.transform(cube)
    return

def myloglike(cube,ndim,nparams):
    if not model.checkSimPriorOK(): return -1.0e99
    realisation=model.realise(cube)
    loglike=poissonLhood(data,realisation)
    return loglike

#-------------------------------------------------------------------------------

def main():

    settingsf='bayestack_settings'
    
    dataset='video'
    survey=surveySetup(dataset)

    binStyle=1
    bins=binSetup(binStyle)

    nlaws=1
    floatNoise=False
    fitter=countModel('ppl',nlaws,settingsf,dataset,binStyle,floatNoise)

    try:
        os.mkdir(outdir)
    except OSError:
        pass

    pymultinest.run(fitter.loglike,fitter.logprior,fitter.nparams,\
                    resume=False,verbose=True,\
                    multimodal=multimodal,max_modes=max_modes,write_output=True,\
                    n_live_points=n_live_points,\
                    evidence_tolerance=evidence_tolerance,\
                    mode_tolerance=-1e90,seed=SEED_SAMP,max_iter=max_iter,\
                    importance_nested_sampling=do_INS,\
                    outputfiles_basename=os.path.join(outdir,outstem),\
                    init_MPI=False)



    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

