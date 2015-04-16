#!/usr/bin/env python

"""
./fluxhisto.py settings.py
"""

import sys
import numpy
import matplotlib.pyplot.as plt

__name_cached=__name__
param_file=sys.argv[-1]

if __name__=='__main__':
    setf='settings'
else:
    setf='%s'%param_file

set_module=importlib.import_module(setf)
globals().update(set_module.__dict__)
__name__=__name_cached


#-------------------------------------------------------------------------------

@profile
def main():
    """
    """

    



#-------------------------------------------------------------------------------

if __name__=='__main__':
    """
    """
    main()
    sys.exit(0)


