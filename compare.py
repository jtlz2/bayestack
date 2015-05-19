#!/usr/bin/env python

"""
This is compare.py
Jonathan Zwart
May 2015

Usage:

./compare.py CHAINS_DIR0 CHAINS_DIR1

"""

import sys
from profile_support import profile
from utils import reportRelativeEvidence

H0=sys.argv[-2]
H1=sys.argv[-1]

#-------------------------------------------------------------------------------


@profile
def main():

    """
    """

    reportRelativeEvidence(H0=H0,H1=H1)


#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
