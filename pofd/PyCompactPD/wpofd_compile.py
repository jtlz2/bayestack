#!/usr/bin/env python

import os,sys

from cffi import FFI
ffibuilder = FFI()

#-------------------------------------------------------------------------------

srcdir='/Users/jtlz2/Dropbox/bayestack/pofd/CompactPD/src'
src='Wpofd.cpp'
src=os.path.join(srcdir,src)
src_str=open(src,'r').read()

os.environ["CC"] = "g++"
libs=['fftw3','m','gsl','gslcblas']# -lfftw3 -lm -lgsl -lgslcblas]
incl=['/usr/local/Cellar/gcc/4.2.1/include/c++/4.2.1/','/usr/local/include/','/usr/local/Cellar/gcc/4.2.1/include/c++/4.2.1/x86_64-apple-darwin13.4.0/',srcdir]

ffibuilder.set_source('_wpofd',src_str,libraries=libs,include_dirs=incl)

ffibuilder.cdef(
"""
struct PD_params{
  double d_max;
  double d_min;
  double source_max;
  double source_min;
  double faint_slop;
  double last_interplot_log10x;
  double last_interplot_log10y;
  double PSFresultionFWHM;
  double pixelsize;
  double sigma_noise;

  double * m_beam;
  int m_beam_size;

  };

double CompactPD_LH(int Nbins, double * DataArray, void * ParamsArray );

""")

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    ffibuilder.compile(verbose=True)
    sys.exit(0)
