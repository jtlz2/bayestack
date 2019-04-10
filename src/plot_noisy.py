#!/usr/bin/env python

"""
Simple simulation demo for Rice distribution

Derived from ipython notebook of the same name

Look for outputs:

rice_counts_QvU.png
rice_counts_input.pdf
rice_counts_recon.pdf

Jonathan Zwart
31 March 2016
"""

import sys
import numpy,scipy
from math import ceil,sqrt
from scipy import integrate
from scipy.interpolate import interp1d
from scipy.stats import rice,rayleigh
import matplotlib.pyplot as plt
from utils import *
from polnUtils import *
from countUtils import writeCountsFile,simulate
from stackUtils import secateur

param_file=sys.argv[-1]
settingsf='%s.bayestack_settings' % param_file

#-------------------------------------------------------------------------------

def main():

    """
    """

    # Import the settings variables
    print 'Settings file is %s' % param_file

    # Import the settings variables
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)

    # Initialize seed
    SEED_SIM=1234
    numpy.random.seed(seed=SEED_SIM)

    # Generate noise-free Q/U data
    # distrib defines dN/dQ and dN/dU, and hence dN/dP0:
    # Q(U) ~ DELTA(Q0(U0))
    # Q(U) ~ U(0,Q0(U0))
    # Q(U) ~ G(Q0(U0),sigma_intrinsic)
    # Q(U) = 0 [OFF]
    N=2981 # was 100000              # Number of objects
    distrib='PPL'           # 'U' or 'G' or 'DELTA' or 'OFF' or 'PPL'
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

    # ... or no signal!
    elif distrib=='OFF':
        Q=numpy.zeros(N)
        U=numpy.zeros(N)

    # ... or a power law signal
    elif distrib=='PPL':
        f='%spost_equal_weights.dat' % outstem
        f=os.path.join(param_file,f)

        # Load equally-weighted posterior samples
        x=numpy.genfromtxt(f)
        nsamp=x.shape[0]
        ncols=x.shape[1]
        summf=os.path.join(param_file,'1-summary.txt')
        summary=numpy.genfromtxt(summf)[-1,:]
        drawmap=summary[-(ncols+1):-2]

        SEED_SIM=7474
#        C=304.0
#        S0=0.12
#        S1=5.8#3.30
#        a0=-1.07#-1.60#-1.07
#        sigma_QU=1.169
        #bins=numpy.linspace(0.0,5.0,51)
   #     C=304.0
   #     S0=0.30
   #     S1=3.75
   #     alpha=-1.58
   #     sigma_QU=1.13#0.8#1.169
#        print BIN_CAT
#        print sigma_QU,S0,S1,C,a0
#        sys.exit(0)
    #    SEED_SIM=7474
    #    C=304.0
    #    S0=0.119
    #    S1=3.30
    #    alpha=-1.08
    #    sigma_QU=1.169
        sigma_QU=S0=S1=S2=S3=C=alpha=beta=gamma=-99.0
        if nlaws==1:
            (sigma_QU,S0,S1,C,alpha)=drawmap
            #(sigma_QU,S0,S1,C,alpha)=(0.106782505590012966E+01,0.224753565298838665E+00,0.295790370249263956E+01,0.462319051197948342E+06,-0.207233157965591541E+01)
            #(sigma_QU,S0,S1,C,alpha)=(0.994914088159975840E+00,0.320033575223691469E+00,0.255445944997367747E+01,0.485991041975387558E+07,-0.182678029396155184E+01)
            #(sigma_QU,S0,S1,C,alpha)=(0.102465519791475335E+01,0.283036175469527918E+00,0.257687943814986742E+01,0.527945739086615845E+03,-0.191352628041741535E+01)
        elif nlaws==2:
            (sigma_QU,S0,S1,S2,C,alpha,beta)=drawmap
        elif nlaws==3:
            (sigma_QU,S0,S1,S2,S3,C,alpha,beta,gamma)=drawmap

        bins=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                        0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0])#,5.8])
        #bins=numpy.linspace(0.0,5.0,26)
        #bins=numpy.linspace(sqrt(2.0)*S0,sqrt(2.0)*S1,101)
        if nlaws!=2:
            pass
            #S1=bins[-1]
       ##     S1*=1.0/sqrt(2.0)
       ##     S0*=1.0/sqrt(2.0)
            #C*=1.0e-10

#        else:
#            S2=bins[-1]
#        print (sigma_QU,S0,S1,S2,S3,C,alpha,beta,gamma)
#        sigma_QU=0.85
        # Parse parameters
        #C=alpha=Smin=Smax=beta=S0=gamma=S1=delta=S2=-99.0
        #values=drawmap
        #nlaws=int(0.5*len(parameters)-1) # infer nlaws from length of param vector
        #params=parameters
        #S1=5.0*SURVEY_NOISE
        print parameters
        print nlaws
        #sys.exit(0)
        if nlaws==1:
            values=[C,S0,S1,alpha]
            params=['C','S0','S1','a0']
        elif nlaws==2:
            values=[C,S0,S1,S2,alpha,beta]
            params=['C','S0','S1','S2','a0','a1']
        elif nlaws==3:
            values=[C,S0,S1,S2,S3,alpha,beta,gamma]
            params=['C','S0','S1','S2','S3','a0','a1','a2']
#        params=['C', 'slope', 'Smin', 'Smax']
#        params=parameters
        #print params
        #print values
        #sys.exit(0)
        #sigma_QU=1.5
        #S1=[values[i] for i in range(len(params)) if params[i]=='S1']
        #iSmax=int([i for i in params if i.startswith('S')][-1][-1])
        #Smax=values[params.index('S%i'%iSmax)]
        nrealns=100
        #sigma_QU=values[params.index('noise')]
        #P0=numpy.zeros(1)
        for i in range(1):
            SEED_SIM+=i
            print values
            print params
            print sigma_QU
            #sys.exit(0)
            print distrib.lower()
            print N,N*nrealns
            #sys.exit(0)
            numpy.random.seed(seed=SEED_SIM)
            Q0=simulate(distrib.lower(),values,params,bins,seed=SEED_SIM,\
                   N=nrealns*N,noise=0.0,dump='Q.txt',output='dummy.txt',\
                   simdocatnoise=False,verbose=True,area=1.0)
            print Q0.min(),Q0.max()
            #sys.exit(0)
            #SEED_SIM+=10
            #SEED_SIM+=i
            #numpy.random.seed(seed=SEED_SIM)
            #if doRayleigh: Q*=0.0 # not necessary?
            #if doRayleigh: sigma_QU=0.110034105777740487E+01
            #sigma_QU=0.000001
            sigma_QU_ray=0.119959401896502138E+01#0.121221464330137763E+01#0.110034105777740487E+01

            nQ=numpy.random.normal(0.0,sigma_QU,nrealns*N)
            nQr=numpy.random.normal(0.0,sigma_QU_ray,nrealns*N)
            Q=Q0+nQ
            #print Q

            U0=simulate(distrib.lower(),values,params,bins,seed=2*SEED_SIM,\
                   N=nrealns*N,noise=0.0,dump='U.txt',output='dummy2.txt',\
                   simdocatnoise=False,verbose=True,area=1.0)
            print U0.min(),U0.max()
            #sys.exit(0)
            #if doRayleigh: U*=0.0 # not necessary?
            nU=numpy.random.normal(0.0,sigma_QU,nrealns*N)
            nUr=numpy.random.normal(0.0,sigma_QU_ray,nrealns*N)
            U=U0+nU
            #print U
            #sys.exit(0)
            P0=numpy.sqrt(numpy.power(Q,2)+numpy.power(U,2))
            Pray=numpy.sqrt(numpy.power(nQr,2)+numpy.power(nUr,2))
            PP=simulate(distrib.lower(),values,params,bins,seed=SEED_SIM,\
                   N=nrealns*N,noise=0.0,dump='PP.txt',output='dummy3.txt',\
                   simdocatnoise=False,verbose=True,area=1.0)
            #PP+=numpy.sqrt(numpy.power(nU,2)+numpy.power(nQ,2))
            PP=rice.rvs(PP/sigma_QU,scale=sigma_QU,size=nrealns*N)
            #print P0
            #numpy.savetxt('P0.txt',P0)
        #P0=rice.rvs(P0/sigma_QU,size=len(P0)) # = OR +=..??
        numpy.savetxt('P0_noisy.txt',P0)
        numpy.savetxt('Pray_noisy.txt',Pray)
        numpy.savetxt('PP_noisy.txt',PP)
        P0=PP
        print sigma_QU

        # Write counts file
        idl_style=False

        # Calculate rescaling factor
        print BIN_CAT
        cat=numpy.genfromtxt(BIN_CAT)

        n=0 # which noise zone?
        #bins=numpy.array([0.03,0.05,0.075,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
        #                              0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,5.0])
        ccat=secateur(cat,cutsDict,n)
        print ccat.shape
        Dcounts=numpy.histogram(Jy2muJy*ccat[:,BIN_COL],bins=bins)[0]
        print Dcounts
        Scounts=numpy.histogram(P0,bins=bins)[0]
        print Scounts
        print bins
        print P0.sum(),Pray.sum()
        #sys.exit(0)
        ScountsRay=numpy.histogram(Pray,bins=bins)[0]
        print ScountsRay,Scounts
        #sys.exit(0)

        #print numpy.argmax(Dcounts)
        #print Dcounts
        #print numpy.max(Dcounts)
        #print Scounts[numpy.argmax(Dcounts)]

        #corrs=numpy.ones(numpy.size(Scounts))*numpy.max(Dcounts)/Scounts[numpy.argmax(Dcounts)]
        print BIN_CAT
        print SURVEY_AREA
        print CORR_BINS
        print bins
   #     sys.exit(0)
        corrs=numpy.ones(Scounts.size)*Dcounts.sum().astype(float)/\
          Scounts.sum().astype(float)
        print corrs
        print CORR_BINS
        print Dcounts.sum(),Scounts.sum()
        #sys.exit(0)
        corrsRay=numpy.ones(ScountsRay.size)*Dcounts.sum().astype(float)/\
          ScountsRay.sum().astype(float)
        #corrs=numpy.ones(numpy.size(Scounts))
        #sys.exit(0)
        print bins
        print P0.min(),P0.max()
        #sys.exit(0)
        print corrs
        print corrsRay

        #sys.exit(0)
        S=writeCountsFile('simP.txt',bins,P0,SURVEY_AREA,\
                          idl_style=idl_style,corrs=corrs)#.astype(float)
        Sray=writeCountsFile('simPray.txt',bins,Pray,SURVEY_AREA,\
                             idl_style=idl_style,corrs=corrsRay)#.astype(float)
        print bins
        print BIN_COL
        print CORR_BINS
        print SURVEY_AREA
        print ccat[:,BIN_COL]
        print BIN_CAT
        print S.size,corrs.size
        print Scounts.sum(),Dcounts.sum()
        #S*=corrs
#        sys.exit(0)
        D=writeCountsFile('datP.txt',bins,Jy2muJy*ccat[:,BIN_COL],SURVEY_AREA,\
            idl_style=idl_style,verbose=False,corrs=CORR_BINS)
        print S
        print Scounts
        print corrs
        print D
        print Dcounts
        print Dcounts.astype(float)/Scounts.astype(float)
        print Dcounts.sum(),Scounts.sum()
#        sys.exit(0)
        # Set up the Rayleigh
        nbins=numpy.size(bins)-1
        ray=numpy.zeros(nbins)
        #mB=medianArray(bins)
        #print sigma_QU_ray,sigma_QU
        #sys.exit(0)
        #sigma_QU_ray=0.121221464330137763E+01#0.110034105777740487E+01#0.119916700779611740E+01#0.120897514423115782E+01#0.123170257903213032E+01#0.123482206906740233E+01#0.120897514423115782E+01#0.143960647406206244E+01
        #if doRayleigh: sigma_QU_ray=sigma_QU
        for ibin in range(nbins):
            #p0=mB[ibin]
            ray[ibin]=rayleigh(scale=sigma_QU_ray).cdf(bins[ibin+1])-\
              rayleigh(scale=sigma_QU_ray).cdf(bins[ibin])
            #ray[ibin]=rices(p0,sigma_QU_ray,bins[ibin],bins[ibin+1],doRayleigh=True)
        print Scounts/Dcounts
        print Dcounts/ray
        print Scounts/ray
#        sys.exit(0)
        corrsRay2=numpy.ones(ray.size)*Dcounts.sum()/ray.sum()
        #corrsRay=numpy.ones(numpy.size(ray))*numpy.max(Scounts*corrs)/\
        #  ray[numpy.argmax(Scounts*corrs)]
        ray*=corrsRay2
        #print ray

        # Now plot the noisy reconstruction
        fig1 = plt.figure(1)
        frame1=fig1.add_axes((.1,.1,.8,.8))
        plt.xlim(0.0,bins[-1]+0.5)#ceil(5.0*sigma_QU))
       # plt.xlim(0.01,10)
        plt.ylim(1.0,600.0)
        if nbins > 50: plt.ylim(1.0,200.0)
        plt.yscale('log')
        #plt.xscale('log')
        plt.xlabel(r'$P_0/\mu$Jy')
        plt.ylabel('Number of objects')
        plt.axvline(1.0*SURVEY_NOISE,color='b',alpha=0.2,label=r'1-5$\sigma$')
        plt.axvline(5.0*SURVEY_NOISE,color='b',alpha=0.2)
        plt.axvline(S0,color='r',alpha=0.2,label=r'MAP range')
        plt.axvline(S1,color='r',alpha=0.2)
        plt.axvline(S2,color='r',alpha=0.2)
        #print S0,S1
        print corrs
        print D.sum(),(S*corrs).sum(),ray.sum()
        plt.errorbar(medianArray(bins),D,fmt='b-',label='noisy data',yerr=numpy.sqrt(D))
        plt.errorbar(medianArray(bins),S*corrs,yerr=numpy.sqrt(S*corrs),fmt='r--',\
                     label='MAP Rice (MC; n%i)'%nlaws)
        plt.errorbar(medianArray(bins),ray,fmt='k--',label='Rayleigh (MC)',yerr=numpy.sqrt(ray))

        plt.errorbar(medianArray(bins),ScountsRay*corrsRay,\
                     yerr=numpy.sqrt(ScountsRay*corrsRay),fmt='g-.',label='Ray (ana)')
        #plt.plot(medianArray(bins),numpy.abs(D-S*corrs),'r.',label='residuals')
        #plt.plot(medianArray(bins),numpy.abs(D-ScountsRay*corrsRay),'g.',label='raysiduals')

        l=plt.legend(loc='upper right',prop={'size':12},frameon=False,numpoints=1)
        linearf='%s/recon_linear_%s.pdf'%(param_file,run_num)

        #Residual plot
        frame2=fig1.add_axes((.1,.1,.8,.1))
        #frame2.set_xticklabels([])
        plt.xlim(0.0,bins[-1]+0.5)
        plt.ylim(-40.0,40.0)
        #plt.yticks(numpy.arange(-50.0,50.0,50.0))
        plt.plot(medianArray(bins),S*corrs-D,'r.-',label='residuals')
        plt.plot(medianArray(bins),ScountsRay*corrsRay-D,'g.-',label='raysiduals')
        #plt.plot(medianArray(bins),ray-D,'k.-',label='ranasiduals')
        plt.axhline(0.0,color='b',alpha=0.2)
        plt.savefig(linearf)
        plt.close()
        print '--> open %s' % linearf
        print
        print 'Model std/sqrtN chisq_red'
        print 'Ray %f %f'%(numpy.std(ScountsRay*corrsRay-D)/numpy.sqrt(D.size-1),calculateReducedChisq(D,ScountsRay*corrsRay,numpy.sqrt(D),ndof=D.size-2))
        print 'Ray2 %f %f'%(numpy.std(ray-D)/numpy.sqrt(D.size-1),calculateReducedChisq(D,ray,numpy.sqrt(D),ndof=D.size-2))
        print 'Rice_n%i %f %f'%(nlaws,numpy.std(S*corrs-D)/numpy.sqrt(D.size-1),calculateReducedChisq(D,S*corrs,numpy.sqrt(D),ndof=D.size-5))
        sys.exit(0)

    # ... or from some other arbitrary distribution we wish to infer....
    ####

    # Now generate noise-free P's
    #Q, U -> P0
    P0=numpy.sqrt(Q**2+U**2)

    # Add gaussian radio noise to each of Q and U,
    # assuming sigma_Q=sigma_U=sigma_QU for now
    sigma_QU=1.164 # was 3.0           #uJy
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

    # Write the P simulation to file
    f='en1jvla/en1jvla_P_a0.txt'
    print 'Minimum flux in catalogue/uJy = %f'%P.min()
    bins=numpy.array([0.02,0.05,0.080,0.10,0.125,0.15,0.20,0.25,0.35,0.5,\
                        0.8,1.0,1.25,1.5,1.75,2.0,2.5,3.0,3.5,4.0,4.5,5.8])
    A=1.979178e-02 # area in sq. deg.
    idl_s=False
    writeCountsFile(f,bins,P,A,idl_style=idl_s,verbose=False)
    sys.exit(0)
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

