#-------------------------------------------------------------------------------

@profile
def calculateDnByDs(bins,counts,eucl=False,verbose=False,idl_style=True,
                    errors=False,bright=False,bright_errors=False,
                    return_all=False):
    """

    This function expects bins to be in uJy
    The output units are uJy^-1, or uJy^1.5 if eucl=True

    The function is area-agnostic: sr^-1 -> sr^-1, deg^-2 -> deg^-2 etc.

    Be careful with units:
       e.g. if bins are in uJy, dn_by_ds will be uJy^-1   = 10^6 Jy
                                dn_ny_ds_eucl in uJy^1.5  = 10^-9 Jy
    etc.
    """

    dn=counts
    Smed=medianArray(bins)
    ds=numpy.absolute(numpy.gradient(Smed))

    if idl_style:
        # i.e. the way DERIV.PRO does it
        ds[0] = abs(-3.0*Smed[0] + 4.0*Smed[1] - Smed[2] / 2.0)
        ds[-1] = abs(3.0*Smed[-1] - 4.0*Smed[-2] + Smed[-3] / 2.0)

    if verbose:
        print 'dn',dn
        print 'ds',ds

    if return_all:
        return dn/ds, (Smed**2.5) * dn/ds,\
          (Smed**2.5) * numpy.sqrt(dn)/ds,\
          (Smed**2.0) * dn/ds,\
          (Smed**2.0) * numpy.sqrt(dn)/ds

    if eucl:
        dn_by_ds = (Smed**2.5) * dn/ds
    elif errors:
        dn_by_ds = (Smed**2.5) * numpy.sqrt(dn)/ds
    elif bright:
        dn_by_ds = (Smed**2.0) * dn/ds
    elif bright_errors:
        dn_by_ds = (Smed**2.0) * numpy.sqrt(dn)/ds
    else:
        dn_by_ds = dn/ds

    return dn_by_ds

#-------------------------------------------------------------------------------


@profile
def powerLawFuncErfsS(S,nlaws,C,alpha,D,beta,Smin,Smax,\
                      Sbinlow,Sbinhigh,S0,gamma,S1,delta,S2,ssigma,area):
    """
    Ketron's equation (5) + (9)

    ; Integration of eq. (7).                                                       
    ;
    ;    /smax               /max_bin                                               
    ;    |     dS * dN/dS *  |        dS' exp(-(S-S')**2/2/rms^2)                      
    ;    /smin               /min_bin                                               
    ;                                                                               
    ; Analytic solution of second integral is error functions.                      
    ;
    """

    if S < Smin or S > Smax:
        return 0.0

    if ssigma > 1.0e-50:
        erfs=0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))
        #erfs *= 0.5
        # Should the next line be min/max or binlow/binhigh???
        # min/max - parameter inferences way off otherwise
        #return erfs * powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Sbinlow,Sbinhigh,S0,area)
        return erfs * powerLawFuncWrap(nlaws,S,C,alpha,D,beta,\
                                       Smin,Smax,S0,gamma,S1,delta,S2,area)
    else: # Noise-free(!!)
        # Should the next line be min/max or binlow/binhigh???
        # THIS NEEDS TO BE THOROUGHLY CHECKED
        return powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Sbinlow,Sbinhigh,\
                         S0,gamma,S1,delta,S2,area)

#-------------------------------------------------------------------------------

@profile
def powerLawFuncWrap(nlaws,S,C,alpha,D,beta,Smin,Smax,S0,gamma,S1,delta,S2,area):

    if nlaws == 1:
        return powerLawFuncS(S,C,alpha,Smin,Smax,area)
    elif nlaws == 2:
        return powerLawFuncDoubleS(S,C,alpha,-99.0,beta,Smin,Smax,S0,area)
    elif nlaws == 3:
        return powerLawFuncTripleS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,area)
    elif nlaws == 4:
        return powerLawFuncQuadS(S,C,alpha,beta,Smin,Smax,S0,gamma,S1,delta,S2,area)

    return

#-------------------------------------------------------------------------------

@profile
def powerLawFuncErrorFn(Si,C,alpha,Smin,Smax,Sbinlow,Sbinhigh,noise,area):

    """
    Units:
            Si, Sbinlow, Sbinhigh, noise - must be the same as each
                                             other, i.e. Jy or uJy
            erfs, alpha - dimensionless
            C - /([Si]*[area]), e.g. Jy^-1 sr^-1, uJy^-1 deg^-2, or some mix
            Hence units of C are tied to units of Si
            area is often just set to 1.0 anyway
    """

    dnds00=True
    if dnds00 and (Si < Smin or Si > Smax):
        #print '***Si',Smin,Smax
        #if False:
        return 0.0

    erfs = erf((Si-Sbinlow)/(sqrt(2.0)*noise)) - erf((Si-Sbinhigh)/(sqrt(2.0)*noise))
    erfs *= 0.5
    n = C * Si**alpha * erfs * area

    return n

#-------------------------------------------------------------------------------

@profile
def calculateI3(C,alpha,Smin,Smax,area,noise=None,dump=None,\
                verbose=False,nlaws=1,D=None,beta=None,S0=None,\
                gamma=-99.0,S1=-99.0,delta=-99.0,S2=-99.0):
    """
    Do this for all bins simultaneously
    Mimic pn_integral.pro:
    ;+                                                                              
    ;                                                                               
    ; Ketron 11/2012                                                                
    ;                                                                               
    ; Integration of eq. (7).                                                       
    ;                                                                               
    ;                                                                               
    ;    /smax               /max_bin                                               
    ;    |     dS * dN/dS *  |        dS' exp(-(S-S')/2/rms^2)                      
    ;    /smin               /min_bin                                               
    ;                                                                               
    ; Analytic solution of second integral is error functions.                      
    ;                                                                               
    ;-  
    """


    #II = numpy.zeros(nbins-1)
    ## Set up the numerical derivative
    #resolution=1000
    #ds = (Smax - Smin) / resolution
    #for ibin in xrange(nbins-1):
    ## This is the integral from Smin to Smax
    ## Step from Smin -> Smax
    #    I=0.0
    #    for istep in xrange(resolution):
    #        Si = Smin + (istep + 0.5) * ds # + -> - because of counting from 0 not 1
    #        erfs = erf((Si-bins[ibin])/(sqrt(2.0)*noise)) - erf((Si-bins[ibin+1])/(sqrt(2.0)*noise))
    #        I += area*C*Si**alpha * erfs * ds
    #    II[ibin]=I

    # I checked this is the same as the numerical integration above
    # (to 1 dec. place.)
    II = numpy.zeros(nbins-1)
    II2 = numpy.zeros(nbins-1)
    # Switch to quadrature integration
    # This is the integral from Smin to Smax
    # Put C scalings here:
    #C *= 10**(1.0/(6.0*alpha))
    #C *= 10**(6.0*alpha)
    #C *= 10**6.0 # Get junk out unless C is of this order (uJy v. Jy)
    for ibin in xrange(nbins-1):
        # This is switched off in order to fix the units of C, but for
        # a long time this was the version in use (1401NNx)
        #II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,area),Smin,Smax)[0]*1000.0

        # I'm trying to investigate the units of the power law in
        # case it affects C
        # **** !!!! This is the wrong place to do it --- the loop compounds!!!!
        # Tried adding a factor of SURVEY_AREA here:
        #print C,alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,area
        #raw_input('pause')
        # 1.0 was = area but I've hard-coded this now
        #C *= 1.0/sqDeg2sr
        #C *= 1.0e6 <- seems to work ok at 140207, but error bars too small
        #C *= 10**(1.0/(6.0*alpha))
        #**** All the quantities are in uJy, except C is in uJy^-1 sr^-1
        #**** alpha is dimensionless
        #**** II needs to be in same area units as ksNoisy (/SURVEY_AREA now)
        sqDeg2srr=sqDeg2sr
        #sqDeg2srr=1.0
        if nlaws == 1:
            # Could 'lambda: S' be deleted - speed-up..?
            # Was Smin,Smax -> Sbinlow,Sbinhigh
            # Need to sort out this maths C expression - II or II2??
            #II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C*10**(1.0/(6.0*alpha)),alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,1.0),bins[ibin],bins[ibin+1])[0]
            # Checked II and II2 are now the same:
            II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C*1.0e-6*10**(-6.0*alpha),alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,sqDeg2srr*SURVEY_AREA),Smin,Smax)[0]
            II2[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,noise/1.0e6,sqDeg2srr*SURVEY_AREA),Smin/1.0e6,Smax/1.0e6)[0]
            #print bins[ibin],bins[ibin+1],Smin,Smax,II[ibin],II2[ibin],II[ibin]/II2[ibin]

            #print C,II[ibin],' '
            #II[ibin]=integrate.quad(lambda S:powerLawFuncErrorFn(S,C,alpha,Smin,Smax,bins[ibin],bins[ibin+1],noise,1.0),Smin,Smax)[0]
            #print II[ibin]

        elif nlaws in [2,3,4]:
            # XXXXX NEED TO FIX INTEGRATION LIMITS XXXXX -> Sbinlow,Sbinhigh
            II[ibin]=integrate.quad(lambda S:powerLawFuncErfsS(S,nlaws,C,alpha,D,beta,Smin/1.0e6,Smax/1.0e6,bins[ibin]/1.0e6,bins[ibin+1]/1.0e6,S0/1.0e6,gamma,S1/1.0e6,delta,S2/1.0e6,noise/1.0e6,sqDeg2srr*SURVEY_AREA),Smin/1.0e6,Smax/1.0e6)[0]

            #print ibin,bins[ibin],bins[ibin+1],Smin,Smax,II[ibin],II2[ibin],II[ibin]/II2[ibin]
        #print II.sum()*1000.0
    #print C,alpha,Smin,Smax
    #sys.exit(0)

    return II

#-------------------------------------------------------------------------------
