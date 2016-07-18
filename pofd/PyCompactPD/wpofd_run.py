#!/usr/bin/env python

import os,sys
import numpy as np
from _wpofd import ffi,lib

VERSION='_v2' # or ''

#-------------------------------------------------------------------------------

def poissonLhood(data,realisation,silent=True):
    if not silent:
        for i in range(len(data)):
            print i,data[i],realisation[i]
    kk=data[np.where(realisation > 0) and np.where(data>0)];
    iii=realisation[np.where(realisation > 0) and np.where(data>0)]
    loglike = (kk*np.log(iii) + kk - kk*np.log(kk) - iii).sum()
    return loglike

#-------------------------------------------------------------------------------

def loglike_pofd(data,model,ParamsArray):
    """
    """

    Nbins=data.shape[0]
    dd=np.ascontiguousarray(data.T)
    DataArray = ffi.cast("double *",dd.ctypes.data)

    xy=model

    i_p=ffi.cast("double *",np.ascontiguousarray(xy).ctypes.data)
    i_l=xy.shape[-1]

    ParamsArray = ffi.new("struct PD_params *",{\
        'd_max':1.0e-3,
        'd_min':10**-7.32,
        'source_max':8.5e-5,
        'source_min':10**-7.32,
        'PSFresultionFWHM':6.0,
        'pixelsize':1.0,
        'sigma_noise':17.0e-6,
        'interplot_length':i_l,
        'interplot_pointer':i_p})

    z=np.zeros(Nbins)
    result=ffi.new("double []",np.ascontiguousarray(z).ctypes.data)
    r=lib.CompactPD_LH(Nbins,DataArray,result,ParamsArray)

    realisation=np.array([result[i] for i in range(Nbins)])

    loglike=poissonLhood(data[:,-1],np.power(10,realisation),silent=True)
    
    return loglike

#-------------------------------------------------------------------------------

def main():

    datadir='/Users/jtlz2/Dropbox/bayestack/pofd/CompactPD%s/Debug/out/'%VERSION
    histof='histo.txt'
    histof=os.path.join(datadir,histof)
    data=np.loadtxt(histof)

    #Npix=data[:,-1].sum()
    #Nbins=data.shape[0]
    #dd=np.ascontiguousarray(data.T)
    #DataArray = ffi.cast("double *",dd.ctypes.data)

    if '2' in VERSION:
        """
        Takes two arrays x, y for the model coordinates
        """
        #x=[-7.32,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9]
        #y=[16.1737,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314]

        # These are very roughly the SKADS counts
        fluxes=[2.0e-8,17.0e-6,85.0e-6]
        dnds_eucl_per_sr=[5.0e-3,1.0,3.0]

        # Some housekeeping
        x=np.log10(np.array(fluxes))
        y=np.log10(np.array(dnds_eucl_per_sr)*np.power(np.power(10,x),-2.5))

        model=np.array([x,y])
        #np.savetxt('xy.txt',model.T)
        i_p=ffi.cast("double *",np.ascontiguousarray(model).ctypes.data)
        i_l=model.shape[-1]

        ParamsArray = ffi.new("struct PD_params *",{\
        'd_max':1.0e-3,
        'd_min':10**-7.32,
        'source_max':8.5e-5,
        'source_min':10**-7.32,
        'PSFresultionFWHM':6.0,
        'pixelsize':1.0,
        'sigma_noise':17.0e-6,
        'interplot_length':i_l,
        'interplot_pointer':i_p})

        #z=np.zeros(Nbins)
        #result=ffi.cast("double *",np.ascontiguousarray(z).ctypes.data)
        #r=lib.CompactPD_LH(Nbins,DataArray,result,ParamsArray)
        #realisation=np.array([result[i] for i in range(Nbins)])
        #loglike=poissonLhood(data[:,-1],np.power(10,realisation),silent=True)
        #print 'loglike = %f [unchecked]\n'%loglike

        loglike=loglike_pofd(data,model,ParamsArray)
        print 'loglike = %f [unchecked]\n'%loglike

    else:
        """
        Takes a data array but uses fixed model (for testing)
        """
        ParamsArray = ffi.new("struct PD_params *",{\
        'd_max':1.0e-3,
        'd_min':np.power(10,-7.32),
        'source_max':8.5e-5,
        'source_min':np.power(10,-7.32),
        'last_interplot_log10x':-7.32,
        'last_interplot_log10y':16.1737,
        'PSFresultionFWHM':6.0,
        'pixelsize':1.0,
        'sigma_noise':17.0e-6})

        result=lib.CompactPD_LH(Nbins,DataArray,ParamsArray)
        print '\nresult: %f [OK]\n' % result
        chisq_true=1553478.027942; tol=1.0e-3
        assert(result-chisq_true<tol), '***FAILURE: chisq is incorrect!!!'

    return 0

#-------------------------------------------------------------------------------

if __name__=='__main__':

    ret=main()
    sys.exit(ret)
