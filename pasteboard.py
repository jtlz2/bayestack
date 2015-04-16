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
