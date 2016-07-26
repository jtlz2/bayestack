#!/usr/bin/env python

import os,sys
from cffi import FFI

from os.path import expanduser
home = expanduser("~")

VERSION='_v2' # or ''

#-------------------------------------------------------------------------------

def main():

    ffibuilder = FFI()

    srcdir='%s/Dropbox/bayestack/pofd/CompactPD%s/src'%(home,VERSION)
    src='Wpofd.cpp'
    src=os.path.join(srcdir,src)
    src_str=open(src,'r').read()

    os.environ["CC"] = "g++"

    libs=['fftw3','m','gsl','gslcblas','stdc++']# -lfftw3 -lm -lgsl -lgslcblas]
    incl=['/usr/local/Cellar/gcc/4.2.1/include/c++/4.2.1/','/usr/local/include/',\
          '/usr/local/Cellar/gcc/4.2.1/include/c++/4.2.1/x86_64-apple-darwin13.4.0/',\
          '/usr/lib/gcc/x86_64-linux-gnu/4.9/include',\
          srcdir]

    ffibuilder.set_source('_wpofd6',src_str,libraries=libs,include_dirs=incl)

    if '2' in VERSION:
        hdr="""
struct PD_params{
  double d_max;
  double d_min;
  double source_max;
  double source_min;
  double faint_slop;
  double PSFresultionFWHM;
  double pixelsize;
  double sigma_noise;

  double * m_beam;
  int m_beam_size;

  int interplot_length;
  double * interplot_pointer;

  };

double CompactPD_LH(int Nbins, double * DataArray, double * result, void * ParamsArray );
"""
    else:
        hdr="""
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

"""
    
    ffibuilder.cdef(hdr)

    ffibuilder.compile(verbose=True)

    return 0

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    ret=main()
    sys.exit(ret)

