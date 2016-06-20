#!/usr/bin/env python

import os,sys
import numpy as np
from _wpofd import ffi,lib

#-------------------------------------------------------------------------------

def main():

    datadir='/Users/jtlz2/Dropbox/bayestack/pofd/CompactPD/Debug/out/'
    histof='histo.txt'
    histof=os.path.join(datadir,histof)
    data=np.loadtxt(histof)

    Nbins=data.shape[0]
    dd=np.ascontiguousarray(data.T)
    DataArray = ffi.cast("double *",dd.ctypes.data)

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
    print '\nresult: %f\n' % result
    chisq_true=1553478.027942; tol=1.0e-3
    assert(result-chisq_true<tol), '***FAILURE: chisq is incorrect!!!'

    return 0

#-------------------------------------------------------------------------------

if __name__=='__main__':

    ret=main()
    sys.exit(ret)
