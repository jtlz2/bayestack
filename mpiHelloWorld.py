#!/usr/bin/env python

"""
Jonathan Zwart
May 2015

Run as

chmod +x mpiHelloWorld.py
mpirun -np NPROCS ./mpiHelloWorld.py
"""

from mpi4py import MPI
import sys
world=MPI.COMM_WORLD

rank=world.rank
size=world.size

print 'MPI processors checked in: rank/size = (%i/%i)\n' % (rank,size)

sys.exit(0)
