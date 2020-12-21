"""
Support functions for bayestack, bayestackClasses

Luminosity functions and evolution thereof

No longer depends on dnds_lumfunc module

Jonathan Zwart
November 2015

"""

import os,sys
import importlib
import glob
import numpy
import mpmath
from math import pi,e,exp,log,log10,isinf,isnan
from scipy import integrate,stats
from scipy.interpolate import interp1d
from scipy.special import erf
import matplotlib.pyplot as plt
from profile_support import profile
from utils import sqDeg2sr,sqrtTwo,find_nearest,medianArray,\
                           interpol,buildCDF,Jy2muJy,interpola
from cosmocalc import cosmocalc

if 'chains' in sys.argv[-1]:
    potential_settings=glob.glob(os.path.join(sys.argv[-1],'*settings*py'))
    assert(len(potential_settings)==1), '***More than one potential settings file!'
    settingsf='.'.join([sys.argv[-1],potential_settings[0].split('/')[-1].split('.')[-2]])
else:
    settingsf=sys.argv[-1].split('.')[-2]

#print '%s is using %s' % (__name__,settingsf)
try:
    set_module=importlib.import_module(settingsf)
    globals().update(set_module.__dict__)
except:
    print '***Warning: Settings not loaded'

#-------------------------------------------------------------------------------
# LF utility functions

@profile
def erfss(S,Sbinlow,Sbinhigh,ssigma):
    """
    For a Gaussian noise model
    """
    if dataset=='sdss':
    	S = float(max(S/1.4,S-0.25))
    else:
	S = float(S)
    return 0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

#-------------------------------------------------------------------------------

def get_Lbins(sbins,z,dl,units='Jy'):
    """
    Units:
        sbins - Jy
        dl - Mpc
        Lbins - W Hz^-1
    """
    if units=='muJy':
       sbins=numpy.array(sbins)*1e-6
    elif units =='mJy':
       sbins=numpy.array(sbins)*1e-3
       
    if isinstance(z,float):
       Lbins = [math.pi*4 *(s*1e-26)* (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1) for s in sbins]
    else:
       Lbins = [math.pi*4 *(s*1e-26)* (dli*3.08e22)**2 * (1 +zi)**((LF_SPEC_INDEX)+1) for s,zi,dli in zip(sbins,z,dl)]
    return numpy.array(Lbins) #W/Hz


#-------------------------------------------------------------------------------

def get_z(ind,z_new=False):
    """
Computes the redshift slice from the binfile

returns z_min, z_max and z_mean
    """
    #z = [ 0.7, 1. , 1.35, 1.7, 2., 2.3, 2.6, 3., 3.5, 4.]
    #z =[0.2, 0.45, 0.7, 1.0, 1.3, 1.6, 1.85, 2.15, 2.35, 2.55, 2.85, 3.15, 3.5]
    if z_new:
        z=[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]
        z_median=[0.22, 0.35, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
    else:   
        z=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]
        z_median=[0.32, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
    #
    #print 'this is num ', ind
    z_min  = z[ind -1]
    z_max  = z[ind]
    print z_min,z_max
    z_mean = numpy.mean((z_min,z_max))
    z_mean = z_median[ind -1]
    
    return z_min,z_max, z_mean
	
def get_num(z_i):	
   '''
   returns the index number corrisponding to z
   '''
   z_i=round(z_i,2)
   z = [0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
   #z_mid =list(numpy.around(medianArray(z),decimals=2))
   z_mid =[0.32, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
   #z_mid=[0.22, 0.35, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
   return z_mid.index(z_i)
#-------------------------------------------------------------------------------

def get_dl(z):
    """
Returns the comoving radial distance in Mpc
    """
    if isinstance(z,float):
       return cosmocalc(z,H0=Ho,WM = wm)['DCMR_Mpc']       
    else:
        
       return [cosmocalc(zi,H0=Ho,WM = wm)['DCMR_Mpc'] for zi in z]
    
#-------------------------------------------------------------------------------

def get_dL(z):
    """
Returns the comoving luminosity distance in Mpc
    """
    if isinstance(z,float):
       return cosmocalc(z,H0=Ho,WM = wm)['DL_Mpc']       
    else:
        
       return [cosmocalc(zi,H0=Ho,WM = wm)['DL_Mpc'] for zi in z]
    
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


    
def get_sbins(fbins,z,dl):
    """
    Units:
        fbins - W Hz^-1 [luminosity bins]
        dl - Mpc
        sbins - ***Jy
    """
    fbins = numpy.array(fbins)
    sbins = fbins/(math.pi*4 * (dl*3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX+1))) #sbins = fbins/(pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    return sbins*1e26 #1e26 #*10.9837 #Jy

#-------------------------------------------------------------------------------
def get_bins(flux_min, flux_max, z, Lbinwidth=0.4, Lmin=20):
        '''
        gives flux bins to use in binner.
        args
        
        flux_min (mJy)  : the minimum (negative) flux in catalogue
        flux_max (mJy)  : Maximum flux in catalogue
        z               : redshift
        Lbin_width(dex) : the binwidth in log10(luminosity) space. it's important because
                          the luminosity function has to be normalized by the binwidth (Mpc^-3 mag^-1).
                          mag corresponds to Lbin_width = 0.4
        Lmin (log10(L)) : the minimum luminosity
        
        returns 
        sbins (mJy)     : the flux bins ready to be used in binner
        '''
        
        dl = get_dl(z)
        Lmax = numpy.log10(get_Lbins([flux_max],z,dl)[0])
        Lbins = numpy.arange(Lmin,round(Lmax+1), Lbinwidth)
        positive_sbins = get_sbins(10**Lbins,z,dl)*1e6
        n = numpy.where(positive_sbins > abs(flux_min))[0][0] +1
        #negative_sbins = 0 - positive_sbins[abs(flux_min) >=positive_sbins]
        negative_sbins =numpy.sort( 0 - positive_sbins[:n])
        sbins = numpy.concatenate([negative_sbins,positive_sbins])
        
        return sbins

#-------------------------------------------------------------------------------

def get_dlds(z,dl):
    """
    dl is comoving distance in Mpc
    dlds is dL/dS in metres^-2
    """
    dlds = (pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
    #print dlds
    return dlds

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

def LFtodnds(lbins, LF, dlds, V, dl,mag):
    """
    Units:
        lbins   - W Hz^-1
        LF(phi) - Mpc^-3 mag^-1
        dlds    - Jy W^-1 Hz
        V(Vmax) - Mpc^3
        dl      - Mpc
    Returns:
        dnds - 
    """

    L = numpy.array(lbins)
    #phi_m = LF*(numpy.log(10**0.4))/L
    phi_m = LF/(numpy.log(10)*L)
    dndl = phi_m * V
    dnds = dndl*dlds
    #print dnds*1e-26
    return dnds*1e-26

def dntoLF(lbins, N, dlds, V):
    """
    Units:
        lbins   - log10 W Hz^-1
        N       - number
        dlds    - Jy W^-1 Hz
        V(Vmax) - Mpc^-3
    Returns:
        LF      - Mpc
    """

    L = numpy.power(10,lbins)
    #nume = N*L/(1e-26*numpy.log(10**0.4))
    nume = N*L*1e26*numpy.log(10**0.4)
    dino = dlds*V
    #print dnds*1e-26
    return nume/dino

#-------------------------------------------------------------------------------

@profile
def schechter(L,ln_Lstar,alpha,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha - no units
        * norm - phi_star Mpc^-3 mag^-1
    Outputs:
        phi
    """
    
    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star*over_Lstar
    phi = norm *(Lr**alpha) *numpy.exp(-Lr)
    #print phi,L, Lstar, alpha, norm,Lr, numpy.exp(-Lr)
    #sys.exit()
    return phi

#-------------------------------------------------------------------------------
def lognormpl(L,ln_Lstar,alpha,sigma,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha - no units
        * norm - phi_star Mpc^-3 mag^-1
        * sigma - no units
    Outputs:
        phi
    """
    
    #print 'parameters'
    
    
    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    sig_1= 1./(2*sigma**2)
    phi = norm *Lr**(1-alpha) * numpy.exp(-sig_1 *numpy.log10(1+ Lr)**2)
    #print phi,L, Lstar, alpha, norm,Lr, numpy.exp(-Lr)
    #print Lstar,alpha, sigma,phi_star
    #sys.exit()
    return phi

@profile
def doublepowerlaw(L,ln_Lstar,alpha1,alpha2,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha1 - no units
        * alpha2 - no units
        * C - phi_0 Mpc^-3 mag^-1
    Outputs:
        phi - Mpc^-3 mag^-1
    """
    
    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    phi = norm* ((Lr)**alpha1 + (Lr)**alpha2)**-1.
    #print L, norm, Lr,phi_star,Lstar, alpha1,alpha2
    #sys.exit()
    return phi

@profile    
def powerlaw(L,ln_Lstar,alpha1,ln_phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha1 - no units
        * C - phi_0 Mpc^-3 mag^-1
    Outputs:
        phi - Mpc^-3 mag^-1
    """

    Lstar    = numpy.power(10,ln_Lstar)
    phi_star = numpy.power(10,ln_phi_star)
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    phi = norm* ((Lr)**alpha1)**-1.
    #print L, norm, Lr,phi_star,Lstar, alpha1,alpha2
    #sys.exit()
    return phi

    
@profile
def doublepower2law(L,Lstar,alpha1,alpha2,phi_star):
    """
    Inputs:
        Lbins - Log10(W Hz^-1)
        * Lstar - W Hz^-1
        * alpha1 - no units
        * alpha2 - no units
        * C - phi_0 Mpc^-3 mag^-1
    Outputs:
        phi - Mpc^-3 mag^-1
    """
    
    over_Lstar = Lstar**-1.
    Lr = L*over_Lstar
    norm=phi_star
    phi = norm* ((Lr)**alpha1 + (Lr)**alpha2)**-1.
    #print L, norm, Lr,phi_star,Lstar, alpha1,alpha2
    #print phi,'does it print once'
    #sys.exit()
    return phi
    
def ddoublepowerlaw(L,ln_phi_star,ln_phi_star2,alpha,ln_L1,beta,ln_L2,gamma,ln_L3,delta):
    '''
    Double double power law
     Inputs:
        Lbins  - Log10(W Hz^-1)
        *L1    - Lstar1,first break [W Hz^-1]
        *L2    - second break 
        *L3    - Lstar2, third break [W Hz^-1]
        *alpha - Slope  
        *beta,gamma,delta - second,third and fouth slope [no units]
        * phi_star - normalization [phi_0 Mpc^-3 mag^-1]
    Outputs:
        phi - Mpc^-3 mag^-1

    '''
    phi_star = numpy.power(10,ln_phi_star)
    phi_star2 = numpy.power(10,ln_phi_star2)
    Lstar1    = numpy.power(10,ln_L1)
    Lstar2    = numpy.power(10,ln_L3)
    L2        = numpy.power(10,ln_L2)  
    
    if L <=L2:
        Lr    = L/Lstar1
        phi   = phi_star* ((Lr)**alpha + (Lr)**beta)**-1.
    elif L2 < L:
        Lr    = L/Lstar2
        phi   = phi_star2* ((Lr)**gamma + (Lr)**delta)**-1.
        
    return phi
#-------------------------------------------------------------------------------
@profile
def LF_SFRD(L, q, params=None, paramsList=None, family=None,z=0):
    """
    Calaclulate SFRD using the LF model
    Cosmology is set globally
    """
    L = 10**L
    sfr = get_sfr_q(1.,q,L)
    if family =='novak':
        L/=(1 + z)**(3.16-z*0.32)
    
    else:
        #L*=(3./1.4)**(-.7)
        #L*=(1.4/3)**(-.7)
        #print 'nani ', family
        L=L
    phi = LF(L=L,params=params,paramsList=paramsList,family=family,SFRD=True,inta=None,area=0,S=0, z=z, dlds=0,Vmax=0,dl=0)
    
    logphi = numpy.log10(phi)+0.
    phi_c = 10**(logphi)
    #print 'disclaimer! This only works for dpl in the faint end'
    if family=='novak' or 'evol' in family:
        return phi*sfr  
        
    return phi_c*sfr


@profile
def dNdS_LF(S, z, dlds, Vmax, dl, mag, params=None, paramsList=None, inta=None, area=None, family=None):
    """
    Source count model
    S is in Jy
    Cosmology is set globally
    """
    L=get_Lbins([S],z,dl)*(1.4/3.)**(-.7)
    area = area*sqDeg2sr
    phi = LF(S, z, dlds,Vmax,dl,params,paramsList,inta,area,family,L[0])
    return LFtodnds(L,phi,dlds,Vmax,dl,mag)

def LF(S, z, dlds, Vmax, dl, params=None, paramsList=None, inta=None, area=None, family=None,L=None,SFRD=False,bypassLim=False,SF=False):

    #Lmin=Lmax=Lnorm=Lstar=Lslope=Lslope2=Lzevol=-99.0

    if family in ['LFsch','LFdpl_dpl','LFdpl_dpl_z','LFdpl_pl','LFpl_dpl','LFlognorm_dpl','LFpl_lognorm']:
        Lnorm=params[paramsList.index('LNORM')]
        Lstar=params[paramsList.index('LSTAR')]
        Lslope=params[paramsList.index('LSLOPE')]
        #Lzevol=params[paramsList.index('LZEVOL')]
    if family in ['LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFdpl_dpl_z','LFlognorm_dpl']:
        Lslope2=params[paramsList.index('LSLOPE2')]
    if family in ['LFlognorm','LFpl','LFpl_dpl','LFdpl','LFdpl_pl', 'LFlognorm_dpl','LFpl_lognorm','LFdpl_dpl','LFdpl_dpl_z']:
        Lnorm_2=params[paramsList.index('LNORM_2')]
        Lstar_2=params[paramsList.index('LSTAR_2')]
        Lslope_2=params[paramsList.index('LSLOPE_2')]

    if family in ['LFlognorm_dpl','LFpl_lognorm','LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFdpl_dpl_z']:
        Lmin2=params[paramsList.index('LMIN2')]
        Lmax2=params[paramsList.index('LMAX2')]
        Lmin2,Lmax2 =numpy.power(10,[Lmin2,Lmax2])

    if family in ['LFdpl_dpl_z','LFevol_dpl','LFevol_logn','LFevol_dpl_L','LFevol_logn_L','LFevol_phi_logn','LFevol_logn_all','LFevol_logn_all_L']:
        alpha_agn=params[paramsList.index('A_agn')]
        alpha_SF=params[paramsList.index('A_SF')]
        beta_agn=params[paramsList.index('B_agn')]
        beta_SF=params[paramsList.index('B_SF')]

    if family in ['LFevol_logn_mat','LFevol_phi_logn_mat','LFevol_logn_el','LFevol_logn_sigma','LFevol_logn_lmin','LFevol_logn_lnorm']:
        alpha_SF=params[paramsList.index('A_SF')]
        alpha_agn=params[paramsList.index('A_agn')]


    if family in ['LFevol_logn_slope','LFevol_logn_all','LFevol_logn_all_L']:
        alpha_SF=params[paramsList.index('A_SF')]
        alpha_agn=params[paramsList.index('A_agn')]
        Lslope2_2=params[paramsList.index('LSLOPE2_2')]

    if family in ['LFevol_logn_lnorm','LFevol_logn_all','LFevol_logn_all_L']:
       	Lnorm=params[paramsList.index('LNORM')]
        Lstar=params[paramsList.index('LSTAR')]
        

    if  'phi' in family:
        alpha_D=params[paramsList.index('A_D')]
        beta_D=params[paramsList.index('B_D')]

    if family in ['LFevol_dpl_s','LFevol_logn_s']:
        alpha_SF=params[paramsList.index('A_SF')]
        beta_SF=params[paramsList.index('B_SF')]

    if family in ['LFdpl', 'LFdpl_dpl','LFdpl_dpl_z']:
        Lslope2_2=params[paramsList.index('LSLOPE2_2')]

    if family in ['LFevol_dpl_L','LFevol_logn_L','LFevol_logn_all_L']:
	   num=get_num(z)
	   #print z, num
	   if num>0:
	       Lmin=params[paramsList.index('LMIN_%d'%num)]
    
    if family in ['LFlognorm','LFpl_lognorm','LFlognorm_dpl','LFevol_logn_sigma','LFevol_logn_all','LFevol_logn_all_L']:
        Lsigma = params[paramsList.index('LSIGMA')]
    #print params
    if L==None:
    	L=get_Lbins([S],z,dl)[0]#*(1.4/3.)**(-0.7)
        
    Lmin=params[paramsList.index('LMIN')]
    Lmax=params[paramsList.index('LMAX')]
    Lmin,Lmax =numpy.power(10,[Lmin,Lmax])
    if bypassLim:
        Lmin=10**15
    #Lmax =numpy.power(10,[Lmax])
    #Lmin = get_Lbins([Lmin/1e6],z,dl)[0]*(1.4/3.)**(-.7)
    #print Lmin,L,Lmax,S,Lnorm_2,Lstar_2,Lslope_2#,Lslope2
    #print Lmin,L,Lmax,S,Lnorm,Lstar,Lslope,Lslope2,Lnorm_2,Lstar_2,Lslope_2,Lsigma
    #print LSTAR_MIN,LSTAR_MAX
    #print family
    #sys.exit()
    if SFRD or SF:
    #This is a little hack to intergate with Lmin and Lmax constraints
        if family in ['novak']:
            return lognormpl(L,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
        elif family in ['LFpl_dpl']:
            return doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
        elif family in ['LFdpl_dpl']:
            return doublepowerlaw(L,Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
        elif family in ['LFpl_lognorm', 'LFlognorm_dpl']:
            #print L,Lstar_2,Lslope_2,Lsigma,Lnorm_2, lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
            return lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
        elif 'evol' in family:
            
            if 'dpl' in family:
                L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                phi_2=doublepowerlaw(L_2,22.41 ,2.4,0.43,-3.50)
                phi= phi_2
            elif 'logn' in family:
                 
                 if 'phi_logn_mat' in family:
                    L_2 = L/(1 + z)**(alpha_SF)
                    phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3))
                    phi_2_evol = (1+z)**(alpha_D)
                    phi = phi_2_evol*phi_2
                    return phi
                    
                 elif 'mat' in family:
                    L_2 = L/(10**(alpha_SF))
                    #print alpha_SF, family
                 elif 'el' in family or '_slope' in family:
                    L_2 = L/(1+z)**(alpha_SF)
                    #phi=lognormpl(L_2,20.8,1.05,0.63, -2.13)
                    #phi_2=lognormpl(L_2,21.36,1.15,0.5, -2.36) #z_6f
                    #phi_2=lognormpl(L_2,21.3,1.14,0.51, -2.32) #z_6g
                    phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3))
                    return phi_2 
                 elif '_lmin' in family:
                    L_2 = L/(1+z)**(alpha_SF)
                    phi=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
                    return phi 
                 elif '_sigma' in family:
                    L_2 = L/(1+z)**(alpha_SF)
                    phi=lognormpl(L_2,numpy.log10(1.85e21),1.22,Lsigma, numpy.log10(3.55e-3))
                    return phi
                 elif '_lnorm' in family:
                    L_2 = L/(1+z)**(alpha_SF)
                    phi=lognormpl(L_2,Lstar,1.22,0.63, Lnorm)
                    return phi
                 elif '_all' in family:
                    L_2 = L/(1+z)**(alpha_SF+z*beta_SF)
                    phi=lognormpl(L_2,Lstar,Lslope2_2,Lsigma, Lnorm)
                    return phi
                 else:
                    L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                    #print alpha_SF, beta_SF, z
                 phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3))
                 #phi_2=lognormpl(L_2,21.16,1.22,0.63, -2.6)
                 phi = phi_2
                 #print L_2 ,numpy.log10(phi), family
                 #print numpy.log10( phi)
            return phi
            
        #print 'no fam ', family
    
    if Lmin < L < Lmax:
        if inta is not None:
            return float(inta(S))
        else:
            #print Smin,Smax,S, Lmin,Lmax, L
            #if family=='LFevol_logn':print 'evol_logn'
            if family=='LFsch':
                #print 'Sch'
                phi=schechter(L,Lstar,Lslope,Lnorm)
               # print 'does it do schecther?'
            elif family=='LFpl':
                phi=powerlaw(L,Lstar_2,Lslope_2,Lnorm_2)            
            elif family=='LFdpl':
                if  Lslope_2 < Lslope2_2:
                #    print numpy.log10(L),'-',           	  
                    return -1.0e99
                phi=doublepowerlaw(L,Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
                #print phi,'does it even come here?????'
            elif family=='LFdpl_pl':
                if  Lslope < Lslope2:
                #print '-',
                    return -1.0e99

                phi_1,phi_2 =0,0
                if Lmin2<L<Lmax:
                    phi_1=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
                if Lmin<L<Lmax2:
                    phi_2=powerlaw(L,Lstar_2,Lslope_2,Lnorm_2)
                phi = phi_1+phi_2

            elif family=='LFpl_dpl':
                if  Lslope < Lslope2:
                #print '-',
                    return -1.0e99

                phi_1,phi_2 =0,0
                lmin2=numpy.log10(Lmin2)
                Lr=10**(lmin2-Lstar)
                Lr_2=10**(lmin2 - Lstar_2)
                Lnorm_2=numpy.log10(10**Lnorm*(Lr_2)**Lslope_2/((Lr)**Lslope + (Lr)**Lslope2))
                #Lnorm_2= Lnorm*(lmin2/Lstar_2)**Lslope_2/((lmin2/Lstar)**Lslope+ (lmin2/Lstar)**Lslope2)
                if Lmin2<L<Lmax:
                    phi_1=powerlaw(L,Lstar_2,Lslope_2,Lnorm_2)
                if Lmin<L<Lmax2:
                    phi_2=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
                phi = phi_1+phi_2



            elif family=='LFdpl_dpl':
		        if Lslope < Lslope2 or Lslope_2 < Lslope2_2:
		        #print '-',
		            return -1.0e99
		        phi_1,phi_2 =0,0
		        if Lmin2<L<Lmax:
		            phi_1=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
		        if Lmin<L<Lmax2:
		            phi_2=doublepowerlaw(L,Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
		        phi = phi_1+phi_2
            
            elif family=='LFdpl_dpl_z':
		        if Lslope < Lslope2 or Lslope_2 < Lslope2_2:
		            return -1.0e99
		        phi_1,phi_2 =0,0
		        if Lmin2<L<Lmax:
		            L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
		            phi_1=doublepowerlaw(L_1,24.29,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
		            #phi_1=doublepowerlaw(L_1,Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
		        if Lmin<L<Lmax2:
		            L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
		            phi_2=doublepowerlaw(L_2,22.11 ,2.4,0.43,-3.10)
		            #phi_2=doublepowerlaw(L_2,,Lstar,Lslope,Lslope2,Lnorm)
		        phi = phi_1+phi_2
            #elif family=='LFevol_logn':print 'evol_logn'    
            elif family in ['LFevol_dpl', 'LFevol_dpl_L', 'LFevol'] :
                    
		            L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
		            L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
		            
		            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
		            phi_2=doublepowerlaw(L_2,22.41 ,2.4,0.43,-3.50)
		            phi = phi_1+phi_2
		            #print S*1e6,z, numpy.log10(L),numpy.log10(L_2),numpy.log10(L_1),numpy.log10(phi_1),numpy.log10(phi_2),numpy.log10(phi)

            elif family=='LFevol_dpl_s':
                        L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                        phi_2=doublepowerlaw(L_2,22.41 ,2.4,0.43,-3.50)
                        phi = phi_2
		            
            elif family in ['LFevol_logn_L']:
		            L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
		            L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)		            
		            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
		            phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
		            phi = phi_1+phi_2
            elif family in [ 'LFevol_logn']:
		            L_1 = L/(1 + z)**(alpha_agn)
		            L_2 = L/(1 + z)**(alpha_SF)		           
		            #L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
		            #L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)		            
		            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
		            phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, -2.6)#numpy.log10(3.55e-3))
		            phi = phi_1+phi_2   
		            

            elif family in ['LFevol_logn_slope'] :
                        L_1 = L/(1 + z)**(alpha_agn)
                        L_2 = L/(1 + z)**(alpha_SF)
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),Lslope2_2,0.63, numpy.log10(3.55e-3))
                        phi = phi_1+ phi_2

            elif family in ['LFevol_logn_sigma'] :
                        L_1 = L/(1 + z)**(alpha_agn)
                        L_2 = L/(1 + z)**(alpha_SF)
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,Lsigma, numpy.log10(3.55e-3))
                        phi = phi_1+ phi_2

            elif family in ['LFevol_logn_lmin'] :
                        L_1 = L/(1 + z)**(alpha_agn)
                        L_2 = L/(1 + z)**(alpha_SF)
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
                        phi = phi_1+ phi_2

            elif family in ['LFevol_logn_lnorm'] :
                        L_1 = L/(1 + z)**(alpha_agn)
                        L_2 = L/(1 + z)**(alpha_SF)
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,Lstar,1.22,0.63, Lnorm)
                        phi = phi_1+ phi_2

            elif family in ['LFevol_logn_all','LFevol_logn_all_L'] :
                        L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
                        L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,Lstar,Lslope2_2,Lsigma, Lnorm)
                        phi = phi_1+ phi_2

   
            elif family=='LFevol_logn_mat':
		            L_1 = L/(10**(alpha_agn))
		            L_2 = L/(10**(alpha_SF))		            
		            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
		            phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
		            #phi_2=lognormpl(L_2,numpy.log10(1.85e21),Lslope2_2,0.63, numpy.log10(3.55e-3))
		            phi = phi_1+phi_2
            elif family=='LFevol_logn_el':
		            L_1 = L/(1+z)**(alpha_agn)
		            L_2 = L/(1+z)**(alpha_SF)		            
		            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
		            #phi_2=lognormpl(L_2,20.8,1.05,0.63, -2.13)
		            #phi_2=lognormpl(L_2,21.36,1.15,0.5, -2.36) #z_6f
		            #phi_2=lognormpl(L_2,21.3,1.14,0.51, -2.32) #z_6g
		            phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3))
		            phi = phi_1+phi_2
            elif family=='LFevol_logn_s':
                        L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
                        phi = phi_2
            elif family in ['LFevol_phi_logn'] :
                        L_1 = L/(1 + z)**(alpha_agn)
                        L_2 = L/(1 + z)**(alpha_SF)
                        #print numpy.log10([L,L_1,L_2]), z
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,20.8,1.05,0.63, -2.13)
                        phi_2_evol = (1+z)**(alpha_D)
                        phi = phi_1+ phi_2_evol*phi_2
            elif family in ['LFevol_phi_logn_mat'] :
                        L_1 = L/(1 + z)**(alpha_agn)
                        L_2 = L/(1 + z)**(alpha_SF)
                        #print numpy.log10([L,L_1,L_2]), z
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3))
                        phi_2_evol = (1+z)**(alpha_D)
                        phi = phi_1+ phi_2_evol*phi_2
      
            elif family=='LFlognorm':
                phi=lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
            elif family=='LFlognorm_dpl':
                #if Lslope < Lslope2:
                    #print '-',
                #    return -1.0e99
                phi_1,phi_2 =0,0
                if Lmin2<L<Lmax:
	                phi_1=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
                if Lmin<L<Lmax2:
                    phi_2=lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
                #print phi_1, phi_2
                phi = phi_1+phi_2
            elif family=='LFpl_lognorm':
                phi_1,phi_2 =0,0
                if Lmin2<L<Lmax:
                    phi_1=powerlaw(L,Lstar,Lslope,Lnorm)
                if Lmin<L<Lmax2:
                    phi_2=lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
                #print phi_1, phi_2
                phi = phi_1+phi_2

            else:print family
            #print ',',L,10**log10phi[0],Lstar,Lslope,Lnorm
            #phi=phi[0]*(1.0+z)**Lzevol
            #if phi==0.0: phi=1.0e-99 # Hmm
            #print phi
            #sys.exit()
            return phi
    else:
    	#print 'passed', Lmin,L,Lmax
        return 0.


#-------------------------------------------------------------------------------

@profile
def IL(dNdS_LF,redshift,dlds,Vmax,dl,mag,params,paramsList,Sbinlow,Sbinhigh,\
       inta=None,area=None,family=None):
    """
    The integrand is the product of dN/dS_LF (count model) and G (the Gaussian)
    Schechter: ['LNORM','LSTAR','LSLOPE','LMIN','LMAX','LZEVOL'] + ['noise']
    Double PL: ['LNORM','LSTAR','LSLOPE','LSLOPE2','LMIN','LMAX','LZEVOL'] + ['noise']
    """
    #print paramsList
    #print params[0],params[1],params[2],params[3]
    Lmin=10**params[paramsList.index('LMIN')]
    Lmax=10**params[paramsList.index('LMAX')]
    [Smin,Smax]=get_sbins([Lmin,Lmax],redshift,dl)
    sigma=params[paramsList.index('noise')]
    #print redshift,dl,Sbinlow,Sbinhigh,Lmin,Lmax,Smin,Smax
    #return integrate.quad(lambda S:dNdS_LF(S,redshift,dlds,Vmax,dl,params=params,paramsList=paramsList,\
    #         inta=inta,area=area,family=family)*erfss(S,Sbinlow/1.0e6,Sbinhigh/1.0e6,sigma/1.0e6),\
    #                                    Smin,Smax,epsabs=1.49e-10, epsrel=1.49e-10,limit=600)[0]

    #intt = mpmath.quad(lambda S:dNdS_LF(S,redshift,dlds,Vmax,dl,mag,params=params,paramsList=paramsList,\
    #         inta=inta,area=area,family=family)[0]*erfss(S,Sbinlow/1.0e6,Sbinhigh/1.0e6,sigma/1.0e6),\
    #                                    [Smin,Smax])
    intt = integrate.quad(lambda x:dNdS_LF(exp(x),redshift,dlds,Vmax,dl,mag,params=params,paramsList=paramsList,\
             inta=inta,area=area,family=family)*erfss(exp(x),Sbinlow/1.0e6,Sbinhigh/1.0e6,sigma/1.0e6)*exp(x),\
                                        log(Smin),log(Smax))[0]
    return intt
#-------------------------------------------------------------------------------

@profile
def calculateL3(params,paramsList,bins=None,area=None,\
                family=None,dump=None,verbose=False,inta=None,dlds = None,redshift=None,Vmax=None,dl=None,mag=None):

    """
    For LF,
    function to calculate mock data for a given power law or interpolation object
    """

    #dl=cosmocalc(redshift,H0=Ho,WM=wm)['DL_Mpc']
    
    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    #print params[1],params[2],params[3],params[4],params[5],params[6]#,params[7]
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        #print bins[ibin],bins[ibin+1]
        #sqDeg2srr=1.0
        #print area,sqDeg2sr
        IIi = abs(IL(dNdS_LF,redshift,dlds,Vmax,dl,mag,params,paramsList,bins[ibin],bins[ibin+1],\
                    inta=inta,area=sqDeg2srr*area,family=family))
        #II[ibin]=abs(IL(dNdS_LF,redshift,dlds,Vmax,dl,params,paramsList,bins[ibin],bins[ibin+1],\
        #            inta=inta,area=sqDeg2srr*area,family=family))
        II[ibin] = IIi
        #print  bins[ibin],bins[ibin+1], IIi

    #print II
    #sys.exit()
    return II
@profile
def get_q(z):
    #Novak 2017
    #z =[0.1, 0.4, 0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]
    #z=numpy.array(z)
    q  =2.78*(1+z)**(-0.14)
    q_l=2.76*(1+z)**(-0.15)
    q_h=2.80*(1+z)**(-0.13)
    
    #mean
    #q=numpy.mean(q)
    #q_l=numpy.mean(q_l)
    #q_h=numpy.mean(q_h) 
    
    #Delhaie et al. 2017
    #q  =2.88*(1+z)**(-0.19)
    #q_l=2.85*(1+z)**(-0.20)
    #q_h=2.91*(1+z)**(-0.18)
    
    #Bell 2013
    #q = 2.64
    #q_l=2.64-0.02
    #q_h=2.64+0.02
    #print 'mean q =',q
    return q_l,q,q_h
    
@profile
def get_sfr_q(F_imf,q,L):
    '''To get the SFR_radio in each bin
       q : Infrared radio correlation
       L : Luminosity range
    '''
 
    SFR = (F_imf*10**-24)*(10**q)*L
    return SFR
    
def get_sfrd_z(phi,logL,sfr, params, paramsList,family,q,l_up,l_low,z=0):

     # phi : the LF in each redshift bin
     #l_low: lower limit of the luminosity range
    #l_up : upper limit of the luminosity range
    #This computes the sfrd in a given redshift bin
   # Apply a spline to get an equation of the line in order to integrate
    #print phi[:30:-40],logL[:30:-40],sfr[0],l_up,l_low
    func = (sfr)*(phi)
    #interp = interp1d(logL,func, kind='cubic')
    #integral = integrate.quad(interp, l_low, l_up)
    #print integral
    intt = integrate.quad(lambda L:LF_SFRD(L,q,params=params,paramsList=paramsList,\
             family=family,z=z),l_low,l_up)
    
       #print intt
    #print intt#, integral
    return intt#integral
    
def sfrd_Behroozi(z,z0=1.243,A=-0.997,B=0.241,C=0.180):
    z1= z - z0
    sfrd = 10**(A*z1) + 10**(B*z1)
    return C/sfrd
    
def sfrd_Madau(z, A=0.015,B=2.7, C=2.9, D=5.6):
    z1 = z + 1.
    top= z1**B
    sfrd = 1. + (z1/C)**D
    return A*top/sfrd/1.65
def parsa_z(z,a=0.18,b=1.04,c=1.77):
    co= a/(b*math.sqrt(2*math.pi))
    ex=-(z-c)**2/(2*b**2)
    return co*numpy.exp(ex)
#-------------------------------------------------------------------------------
