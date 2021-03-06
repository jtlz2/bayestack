#!/usr/bin/env python

'''
This is dnds_lumfunc.
Eliab Malefahlo
October 2015
 
This function converts dnds into LF
Usage:
    ./dnds_lumfunc.py CHAINS_DIR

Arguments
dnds  : The reconstructed source counts (from reconstruct)
sbins : The flux(uJy) bins that coresponds to dnds
z_min : minimum redshift of slice
z_max : max redshift of slice

returns 
rho_m : the normalized luminosity function (Log10[Mpc^-3 mag^-1])
Lbins : Log Luminosity bins corresponding to rho_m (Log10[L])

dnds = [4.78118760e-05,   6.33563041e-05,   8.00873144e-05,   8.50948383e-05]
sb = [ 45.,     55. ,    65.,     75.]
dnds_lumfunc.get_lumfunc(dnds,sb,1.8,2.5)

output
(array([ -7.60008913e+00,  -7.52893349e+00,  -7.46877069e+00,
         1.21271141e-33]),
 array([ 24.35010723,  24.43725741,  24.50980808,  24.57195598]))

'''

#from numpy import *
import numpy
import os,sys,math,shutil
import importlib
from pylab import*
from matplotlib.ticker import AutoMinorLocator
from cosmocalc import cosmocalc
#import pymultinest
from bayestackClasses import countModel
from utils import sqDeg2sr,fetchStats
from countUtils import calculateDnByDs,medianArray

param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file
Ho=71.0
wm=0.27
LF_SPEC_INDEX = -0.7

#expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
#print expt.data
#print len(expt.bins)
#print expt.bins
# Fetch best-fit parameters and calculate best-fit line
#plotTruth=dict((name,-99.0) for name in expt.parameters)
#stats=fetchStats(outdir,expt.parameters,plotTruth)
#SMIN_MAP=stats['S0'][0]
#SMIN_MAP_UPPER= stats['S0'][-2]
#SMIN_MAP_LOWER=stats['S0'][-1]



#-------------------------------------------------------------------------------

def get_Vmax(zlo,zup):
    z  = zup
    V1 = cosmocalc(z,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
    V2 = cosmocalc(zlo,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
    area = SURVEY_AREA*sqDeg2sr
    vmax = area*(V1-V2)/(4*pi)
    #print 'vmax',vmax,V1,V2,zup,zlo,SURVEY_AREA*sqDeg2sr,sqDeg2sr
    return vmax

#-------------------------------------------------------------------------------

    
def get_dnds(dnds_ecl,sbins):
    """
    Maps dnds_eucl -> dnds
    Units:
        sbins - muJy
        dnds_eucl - standard units
    """
    dnds=[]
    for i in range(len(sbins)):
        s2_5 = 1e-26*(sbins[i]*1e-6)**2.5
        print s2_5,dnds_ecl[i],SURVEY_AREA*sqDeg2sr
        dndsi=(dnds_ecl[i]*SURVEY_AREA*sqDeg2sr)/(s2_5)
        dnds.append(dndsi)
    return numpy.array(dnds)

#-------------------------------------------------------------------------------


def get_dnds_ecl(dnds,sbins):
    """
    Maps dnds -> dnds_eucl
    Units:
        sbins - muJy
        dnds - standard units
    """
    s2_5 = 1e-26 *(sbins)**2.5
    print 's^2.5'
    print s2_5
    dnds_ecl = dnds *s2_5/SURVEY_AREA/sqDeg2sr
    return dnds_ecl

#-------------------------------------------------------------------------------

def get_z(ind):
    """
Computes the redshift slice from the binfile

returns z_min, z_max and z_mean
    """
    #z = [ 0.7, 1. , 1.35, 1.7, 2., 2.3, 2.6, 3., 3.5, 4.]
    z =[0.2, 0.45, 0.7, 1.0, 1.3, 1.6, 1.85, 2.15, 2.35, 2.55, 2.85, 3.15, 3.5]
    
    z_min  = z[ind -1]
    z_max  = z[ind]
    z_mean = mean((z_min,z_max))
    
    return z_min,z_max, z_mean
    #return 1.8,2.5,2.15
	
#-------------------------------------------------------------------------------

def get_dl(z):
    """
Returns the luminosity distance in Mpc
    """
    return cosmocalc(z,H0=Ho,WM = wm)['DL_Mpc']
    
#-------------------------------------------------------------------------------

def get_dsdl(z,dl):
    """
    dl is luminosity distance in Mpc
    dsdl is dS/dL in metres^-2
    """
    dsdl = 1./(math.pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    print 'dsdl'
    print dsdl
    return dsdl

#-------------------------------------------------------------------------------


def get_Lbins(sbins,z,dl):
    """
    Units:
        sbins - muJy
        dl - Mpc
        Lbins - W Hz^-1
    """
    Lbins = [math.pi*4 *(s*1e-32)* (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1) for s in sbins]
    return Lbins #W/Hz

#-------------------------------------------------------------------------------

    
def get_sbins(fbins,z,dl):
    """
    Units:
        fbins - W Hz^-1 [luminosity bins]
        dl - Mpc
        sbins - ***Jy
    """
    sbins = fbins/(math.pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    return sbins*1e26 #Jy

#-------------------------------------------------------------------------------

    
def LFtodnds(lbins, LF , z_min, z_max):
    """
    Units:
        lbins - Log10(W Hz^-1)
        LF - rho_m Mpc^-3 mag^-1
    Returns:
        sbins - flux bins in muJy
        dnds_ecl - dnds_eucl in standard units

    Converts LF(rho_m) to source coundts (dn/ds s^2.5)
    returns sbins, dnds
    """

#    Lbins = medianArray(lbins)
    LF=numpy.power(10,LF)
    Lbins = numpy.power(10,lbins)
    V = get_Vmax(z_min,z_max) 
    print 'Vmax'
    print V
    mag = [(lbins[i+1] - lbins[i])/0.4 for i in range(len(lbins)-1)]
    mag.append(mag[-1])
    
    print 'mag'
    print mag
    z    = mean((z_min,z_max))
    dl   = get_dl(z)
    dsdl = get_dsdl(z,dl)
    dndl = LF * V*mag/Lbins
    print 'dndl'
    print dndl
    sbins = get_sbins(Lbins,z,dl)
    print 'sbins'
    print sbins
    dnds = dndl/dsdl
    print 'dnds'
    print dnds
    dnds_ecl = get_dnds_ecl(dnds,sbins)
    
    return sbins*1e6,dnds_ecl

#-------------------------------------------------------------------------------

        
def get_lumfunc(dnds_ecl,sbins,z_min,z_max,dnds_ecl_up=[],dnds_ecl_down=[]):
    """
    Units:
        dnds_ecl 		- dnds_eucl in standard units
        sbins 			- flux bins in muJy
        dnds_ecl_up 	- upper error bound
        dnds_ecl_down 	- lowwer error bound 
    Returns:
        Lbins 		- Log10(W Hz^-1)
        Log_rho_m (LF) 	- log(rho_m Mpc^-3 mag^-1)
        log_err_up		- upper bound of LF 
        log_err_down	- lower bound of LF 
        
    Converts source counts (dnds s^2.5) to LF (rho_m)
    returns Lbins, LF 
    """

    rho_m      = numpy.zeros(len(sbins)-1)
    dnds_ecl = numpy.array(dnds_ecl)
    sbins    = medianArray(sbins)
    print numpy.where(sbins>0)
    
    dnds_ecl = dnds_ecl[numpy.where(sbins>0)]
    sbins    = sbins[numpy.where(sbins>0)]
    dnds_ecl_up = dnds_ecl_up[numpy.where(sbins>0)]
    dnds_ecl_down =dnds_ecl_down[numpy.where(sbins>0)]
    print 'recon'
    print dnds_ecl
    
    #Luminosity distance
    z  = mean((z_min,z_max))
    dl = cosmocalc(z,H0=Ho,WM = wm)['DL_Mpc'] #mean LDist
    dl_min =cosmocalc(z_min,H0=Ho,WM = wm)['DL_Mpc'] #minimum Ldist 
    dl_max =cosmocalc(z_max,H0=Ho,WM = wm)['DL_Mpc'] #maximu
    
    #error esitmate
    dl_err_down = dl - dl_min
    dl_err_up   = dl_max - dl
    err_up = numpy.sqrt( (dl_err_up/dl)**2 + (dnds_ecl_up/dnds_ecl)**2 )
    err_down = numpy.sqrt( (dl_err_down/dl)**2 + (dnds_ecl_down/dnds_ecl)**2 )
    print 'errors\n'
    print dl_err_down, dl_err_up, err_up, err_down,'\n'
    Vmax=get_Vmax(z_min,z_max)
    o_Vmax    =  1./(Vmax)
    print '1/vmax'
    print o_Vmax,SURVEY_AREA*sqDeg2sr
    dnds = get_dnds(dnds_ecl,sbins)
    
    print 'dsdn'
    print dnds
    dndl = get_dsdl(z,dl)*dnds
    rho     = dndl*o_Vmax
    
    print 'dndl'
    print dndl
    print 'rho'
    print rho_m

    L = get_Lbins(sbins,z,dl) 
    Lbins =log10(L)
    
    dL   = [(Lbins[i+1] - Lbins[i])/0.4 for i in range(len(Lbins)-1)]
    dL.append(dL[-1]) #adding last dL for dimensions
    
    rho_m     = rho*L/dL
    log_rho_m = log10(rho_m)
    rho_err_up = rho_m*err_up
    rho_err_down = rho_m*err_down
		
    log_err_up = numpy.log10(rho_m + rho_err_up) - numpy.log10(rho_m)
    log_err_down = numpy.log10(rho_m) - numpy.log10(rho_m - rho_err_down)  
	
    print 'more errors \n \n'
    print rho_err_up, rho_err_down,rho_m ,log_err_up, log_err_down,'\n \n'

	
    print log_rho_m,Lbins
    if len(dnds_ecl_up) > 0.:
    	return Lbins,log_rho_m, log_err_up, log_err_down
    return Lbins,log_rho_m

#-------------------------------------------------------------------------------


def schechter(Lbins,Lstar,alpha, norm):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        Lstar - W Hz^-1
        alpha - no units
        norm - phi_star Mpc^-3 mag^-1
    Outputs:
        Lbins -
        log10(phi) -
    """
    if type(Lbins) != float:
        Lbins = medianArray(Lbins)

    Lbins = numpy.power(10,Lbins)
    Lr = Lbins/Lstar
    print Lbins, Lstar, alpha, norm
    phi = norm *(Lr)**alpha *numpy.exp(-Lr)/Lstar
    return Lbins, log10(phi)

#-------------------------------------------------------------------------------
def doublepowerlaw(Lbins,Lstar,alpha1, alpha2, C):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        Lstar - W Hz^-1
        alpha - no units
        norm - phi_star Mpc^-3 mag^-1
    Outputs:
        Lbins -
        log10(phi) -
    """
    if type(Lbins) != float:
        Lbins = medianArray(Lbins)

    Lbins = numpy.power(10,Lbins)
    Lr = Lbins/Lstar
    print Lbins, Lstar, alpha1,alpha2, C
    phi = C/ ((Lr)**alpha1 + (Lr)**alpha2)
    return Lbins, log10(phi)
    
#-------------------------------------------------------------------------------

def LFToDnByDs(Lbins,LFtype,LFparams,zlow,zhigh):
    """
    Wrapper function to convert lumfunc parameters to source counts
    """
    if LFtype=='LFschechter':
        Lstar,alpha,norm=LFparams
        Lbins,phi=schechter(Lbins,Lstar,alpha,norm)
        Sbins,dNdS=LFtodnds(Lbins,phi,zlow,zhigh)
    elif LFtype=='LFdoublePL':
        pass

    return Sbins,dNdS

#-------------------------------------------------------------------------------

def testLumFuncs():


    s=loadtxt('%s/recon_stats.txt'%outdir)

    z = [ 0.7, 1. , 1.35, 1.7, 2., 2.3, 2.6, 3., 3.5, 4.]
    #determining the redshift from the name of the datafile
    if 'sdss' in datafile:
        ind = int(datafile[-5])
        z_min= z[ind -1]
        z_max= z[ind]
        print z_min,z_max

    # Testing get_lumfunc and LFtodnds
    dnds = [4.78118760e-05,   6.33563041e-05,   8.00873144e-05,   8.50948383e-05]
    sb = [ 45.,     55. ,    65.,     75.]
    lbins,LF=get_lumfunc(dnds,sb,z_min,z_max)
    sb_out,dnds_out=LFtodnds(lbins,LF,z_min,z_max)
    s=loadtxt('%s/recon_stats.txt'%outdir)
    xrecon=s[:-1,0]; yrecon=s[:-1,1]
    yrecon_down=s[:-1,2]; yrecon_up=s[:-1,3]
            
    L,rho     = get_lumfunc(yrecon,xrecon,z_min,z_max)
    
########## parameters##########    
    phi =5e16
    Ls = 6e25
    alpha=-1
    a='alpha'
    alpha1 = 3.
    alpha2 = 1.3
    C = 8e-8
#_____________________________#    
    
    L_s,rho_s = schechter(L,Ls,alpha, phi)
    L2,rho2 = doublepowerlaw(L,Ls,alpha1,alpha2,C)
    
    plot(L,rho,'s',label='%3.2f < z < %3.2f'%(z_min,z_max))
    plot(L[:-1],rho_s,label='Schechter')
    #plot(L[:-1],rho2,label='Double power-law')
    tick_params(axis='both',which = 'major', labelsize=15,width =2)
    tick_params(axis='both',which = 'minor', labelsize=12, width=1)
    text(23,-8.9,'  $\phi_* =%5.2e$ \n  $L_* =%5.2e$ \n  $\%s =%5.2f$'%(phi,Ls,a,alpha),fontsize=17)
    xlabel('1.4 GHz log[L(W/Hz)]',fontsize = 18)
    ylabel(r'log[$\rho_m$(mpc$^{-3}$mag$^{-1}$)]',fontsize = 18)

    
    #ax.xaxis.set_minor_locator(AutoMinorLocator())
    #ax.yaxis.set_minor_locator(AutoMinorLocator())
    
    legend().draggable()
    show()
	
    #xscale('log')
    #yscale('log')
    #show()
    lff='lf.pdf'
    savefig(lff)
    print '--> Look in %s' % lff	

    assert(numpy.allclose(dnds_out,dnds[:-1])), '**get_lumfunc <-> LFtodnds do not match'
    print '***All tests passed OK'

    # Test Schechter
    
    return

#-------------------------------------------------------------------------------

def main():
    # Import the settings variables
    print 'Settings file is %s' % param_file

    # Import the settings variables
    set_module=importlib.import_module(setf)
    globals().update(set_module.__dict__)

    testLumFuncs()
    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)

#-------------------------------------------------------------------------------
