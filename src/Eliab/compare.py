#!/usr/bin/env python

"""
This is compare.py
Jonathan Zwart
May 2015

Usage:

./compare.py chains_150520?
./compare.py chains_15*
./compare.py chains_150520a chains_150520b ... chains_150520c
etc.

"""

import sys
from profile_support import profile
from utils import reportRelativeEvidences, myRelativeEvidences

globList=sys.argv[1:]

#-------------------------------------------------------------------------------


@profile
def main():

    """
    """

    #reportRelativeEvidence(H0=H0,H1=H1)
    assert(len(globList)>1), '***Require >1 hypothesis for model comparison!!'
    reportRelativeEvidences(globList)


#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
