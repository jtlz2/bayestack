#!/bin/bash
#PBS -N 160722gbl
#PBS -o $HOME/Dropbox/bayestack/pofd/PyCompactPD/chains_test6/output-$PBS_JOBID.log
#PBS -e $HOME/Dropbox/bayestack/pofd/PyCompactPD/chains_test6/output-$PBS_JOBID.err
#PBS -S /bin/bash
#PBS -l walltime=72:00:00
#PBS -l nodes=1:ppn=48
#PBS -l mem=100MB
##PBS -m bea
##PBS -M jtlz2@astro.columbia.edu

ulimit
#ulimit -u 400000

export OMP_NUM_THREADS=1

RUNPATH=$HOME/Dropbox/bayestack/pofd/PyCompactPD
cd $RUNPATH
date
pwd

./wpofd_compile.py
date
cd $RUNPATH
pwd
#LD_LIBRARY_PATH="/usr/local/lib/" mpiexec.osc -n $PBS_NP ./wpofd_fit.py
LD_LIBRARY_PATH="/usr/local/lib/" valgrind --tool=massif --time-unit=B --massif-out-file=./mem_jz.out ./wfit.sh
#LD_LIBRARY_PATH="/usr/local/lib/" valgrind ./wfit.sh

date

