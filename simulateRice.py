#!/usr/bin/env python

"""
Simple simulation demo for Rice distribution

Derived from ipython notebook of the same name

Look for outputs:

rice_counts_QvU.png
rice_counts_input.pdf
rice_counts_recon.pdf

Jonathan Zwart
16 September 2015
"""

import sys
import numpy,scipy
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.stats import rice,rayleigh
import matplotlib.pyplot as plt
from utils import *
from polnUtils import *

#-------------------------------------------------------------------------------

def main():

    # Initialize seed
    SEED_SIM=1234
    numpy.random.seed(seed=SEED_SIM)

    # Generate noise-free Q/U data
    # distrib defines dN/dQ and dN/dU, and hence dN/dP0:
    # Q(U) ~ DELTA(Q0(U0))
    # Q(U) ~ U(0,Q0(U0))
    # Q(U) ~ G(Q0(U0),sigma_intrinsic)
    N=100000              # Number of objects
    distrib='U'           # 'U' or 'G' or 'DELTA'
    Q0=20.0; U0=10.0      # uJy
    sigma_intrinsic=10.0  # uJy ***NB THIS IS **NOT** THE RADIO NOISE

    # ***For now, assume all objects have same underlying Q,U
    if distrib=='DELTA':
        Q=Q0*numpy.ones(N)
        U=U0*numpy.ones(N)

    # ... or are from a uniform distribution
    elif distrib=='U':
        Q=Q0*numpy.random.rand(N)
        U=U0*numpy.random.rand(N)

    # ... or a Gaussian (peak=0)
    elif distrib=='G':
        Q=numpy.random.normal(0.0,sigma_intrinsic,N)
        U=numpy.random.normal(0.0,sigma_intrinsic,N)

    # ... or from some other arbitrary distribution we wish to infer....
    ####

    # Now generate noise-free P's
    #Q, U -> P0
    P0=numpy.sqrt(Q**2+U**2)

    # Add gaussian radio noise to each of Q and U,
    # assuming sigma_Q=sigma_U=sigma_QU for now
    sigma_QU=3.0           #uJy
    Q+=numpy.random.normal(0.0,sigma_QU,N)
    U+=numpy.random.normal(0.0,sigma_QU,N)

    # Plot Q-U plane
    plt.xlim(-3.0*Q0,3.0*Q0)
    plt.ylim(-3.0*Q0,3.0*Q0)
    plt.xlabel('Q')
    plt.ylabel('U')
    plt.plot(Q,U,'.',label='$Q$ v. $U$')
    plt.plot(Q0,U0,'r.',label='($Q_0$,$U_0$)')
    plt.plot(Q,Q,'k-',label='$Q=U$')
    l=plt.legend(frameon=False,numpoints=1)
    plt.savefig('rice_counts_QvU.png')
    plt.close()

    # Generate (noisy) mock P data (given noisy Q and U only, i.e. as in real life)
    # noisy Q,U -> noisy P
    P=numpy.sqrt(Q**2+U**2)

    # From the noise-free P0, add Rician noise
    # NB http://www.mathworks.com/matlabcentral/fileexchange/14237-rice-rician-distribution/content/rician/html/ricedemo.html
    #Pprime=rice.rvs(P,P0,scale=sigma_QU,size=N)
    # P0 -> P
    Pprime=rice.rvs(P0/sigma_QU,scale=sigma_QU,size=N) # (*)

    # This plot is f(Q0,U0) [U, DELTA, G, etc.] + n -> f(P0) + n ->
    # g(P)|f(P0), showing:
    # 1. Underlying noisy Q/U distributions
    # 2. Noisy measured P distribution (generated from (1))
    # 3. Rician Pprime corresponding to the noise-free data
    #    (i.e. a function of P0 and its distribution dN/dP0)
    # NB The plot shows (3) is a fair tracer of (2), the noisy data, so we can use (*) above to
    #    "add rician noise" to any noise-free distribution dN/dP0,
    #    at different S/N and for different distributions of Q0/U0

    # These bins define the plotting and the stored noise-free distribution dN/dP0
    bins=numpy.linspace(0.0,2.0*max(Q0,U0),101)
    plt.xlabel(r'S/$\mu$Jy')
    plt.ylabel('number of objects')
    plt.hist(Q,bins=bins,color='k',label=r'Q=%3.1f$\pm$%3.1f'%(Q0,sigma_QU),alpha=0.1)
    plt.hist(U,bins=bins,color='g',label=r'U=%3.1f$\pm$%3.1f'%(U0,sigma_QU),alpha=0.1)
    plt.hist(P,bins=bins,color='b',alpha=0.5,\
             label=r'P$\equiv\sqrt{Q^2+U^2}$') # noisy data
    plt.hist(Pprime,bins=bins,color='r',alpha=0.5,\
             label=r'Rice(%3.1f,%3.1f)'%(numpy.sqrt(Q0**2+U0**2),sigma_QU))
    title='%s (Q0=%3.1f, U0=%3.1f, sig=%3.1f, N=%i)'%(distrib,Q0,U0,sigma_intrinsic,N)
    plt.title(title)

    # n0 and b0 specify dN/dP0
    n0,b0,pp=plt.hist(P0,bins=bins,color='g',label='P0 noise-free')
    legend=plt.legend(frameon=False)
    plt.savefig('rice_counts_input.pdf')
    plt.close()

    # Test some integrals
    # http://docs.scipy.org/doc/scipy/reference/tutorial/integrate.html#general-multiple-integration-dblquad-tplquad-nquad
    area = integrate.dblquad(lambda x, y: x*y, 0, 0.5, lambda x: 0, lambda x: 1-2*x)
    print area
    area2 = integrate.dblquad(lambda
                              p,p0:1.0*F(p,p0,1.0*sigma_QU),1.0,20.0,\
                              lambda p:2.0,lambda p:2.5)
    print area2

    # Now we want to generate a function that mimics the simulated
    # histogram for noisy P above
    sigma_QU=3.0 # uJy Noise on Q/U

    # Evalulate the function at these points:
    nbins=10
    binslo=0.0; binshi=30.0 # uJy
    Pbins=numpy.linspace(binslo,binshi,nbins+1)

    # Model parameters
    # dN/dP0 distribution object
    inta=P0Dist(b0,n0)
    #inta=None # Toggle this to enter powerLaw mode
    # Start and stop range
    Pmin=1.0e-4 # uJy
    Pmax=40.0 # or e.g. numpy.sqrt(Q0**2+U0**2)

    paramsList=['S0','S1','noise','C','a0']
    params=[Pmin,Pmax,sigma_QU,1.0,-1.0]

    # Calculate the noisy integral
    II2=calculateP3(params,paramsList,bins=Pbins,area=1.0/sqDeg2sr,\
                    family=None,dump=None,verbose=False,inta=inta)

    plt.xlim(binslo,binshi+5.0)
    #plt.yscale('log')
    #plt.plot(medianArray(Pbins),numpy.vectorize(dNdP0)(medianArray(Pbins),Pmin,Pmax),'r')
    plt.xlabel(r'S/$\mu$Jy')
    plt.ylabel('number of objects')
    plt.axvline(Pmin,alpha=0.25)
    plt.axvline(Pmax,alpha=0.25)
    n0,b0,pp=plt.hist(P0,bins=bins,color='g',label='P0 noise-free')
    u,v,w=plt.hist(Pprime,bins=bins,color='r',label='Pprime noisy')
    x,y,z=plt.hist(P,bins=bins,color='b',alpha=0.5,label=r'P$\equiv\sqrt{Q^2+U^2}$')

    rescale=x.max()/II2.max()
    plt.plot(medianArray(Pbins),rescale*II2,'k',label='II2 (analytic)')

    l=plt.legend(frameon=False)
    plt.savefig('rice_counts_recon.pdf')
    plt.close()

    print 'Finished!!!'

    return 0

#-------------------------------------------------------------------------------

if __name__=='__main__':
    ret=main()
    sys.exit(ret)

#-------------------------------------------------------------------------------

