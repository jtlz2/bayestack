#!/usr/bin/env python

import os,sys
import numpy as np
from _wpofd import ffi,lib

VERSION='_v2' # or ''

#-------------------------------------------------------------------------------

def main():

    datadir='/Users/jtlz2/Dropbox/bayestack/pofd/CompactPD%s/Debug/out/'%VERSION
    histof='histo.txt'
    histof=os.path.join(datadir,histof)
    data=np.loadtxt(histof)

    Nbins=data.shape[0]
    dd=np.ascontiguousarray(data.T)
    DataArray = ffi.cast("double *",dd.ctypes.data)
    #xprint dd
    if '2' in VERSION:
        """
        Takes two arrays x, y for the model coordinates
        """
        x=[-7.32,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9]
        y=[16.1737,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314]
        #y=range(10,18)

        xy=np.array([x,y])
        i_p=ffi.cast("double *",np.ascontiguousarray(xy).ctypes.data)

        ParamsArray = ffi.new("struct PD_params *",{\
        'd_max':1.0e-3,
        'd_min':10**-7.32,
        'source_max':8.5e-5,
        'source_min':10**-7.32,
        'PSFresultionFWHM':6.0,
        'pixelsize':1.0,
        'sigma_noise':17.0e-6,
        'interplot_length':8,
        'interplot_pointer':i_p})

        z=np.zeros(Nbins)
        result=ffi.cast("double *",np.ascontiguousarray(z).ctypes.data)
        r=lib.CompactPD_LH(Nbins,DataArray,result,ParamsArray)

        loglike=result[0] # Dereferencing trick, i.e. *p
        print '\nresult: %f [unchecked]\n' % loglike

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
