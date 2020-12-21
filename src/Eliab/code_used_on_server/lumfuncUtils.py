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
    #potential_settings=glob.glob(os.path.join(sys.argv[-1][:sys.argv[-1].index('/')],'*settings*py'))
    potential_settings=glob.glob(os.path.join(sys.argv[-1],'*settings*py'))
    #print 'setting file ',potential_settings[0][:-3]
    assert(len(potential_settings)==1), '***More than one potential settings file!'
    #settingsf=potential_settings[0][:-3]
    settingsf='.'.join([sys.argv[-1],potential_settings[0].split('/')[-1].split('.')[-2]])
else:
    settingsf=sys.argv[-1].split('.')[-2]

#print '%s is using %s' % (__name__,settingsf)
try:
    set_module=importlib.import_module(settingsf)
    print settingsf
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
    	S = float(max(S/1.4,S-2.5e-4))
    else:
	S = float(S)
    return 0.5*(erf((S-Sbinlow)/(sqrtTwo*ssigma)) - erf((S-Sbinhigh)/(sqrtTwo*ssigma)))

#-------------------------------------------------------------------------------

def get_Lbins(sbins,z,dl):
    """
    Units:
        sbins - Jy
        dl - Mpc
        Lbins - W Hz^-1
    """
    if isinstance(z,float):
       Lbins = [math.pi*4 *(s*1e-26)* (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1) for s in sbins]
    else:
        Lbins = [math.pi*4 *(s*1e-26)* (dli* 3.08e22)**2 * (1 + zi)**((LF_SPEC_INDEX)+1) for s,zi,dli in zip(sbins,z,dl)]
    return Lbins #W/Hz


#-------------------------------------------------------------------------------

def get_z(ind):
    """
Computes the redshift slice from the binfile

returns z_min, z_max and z_mean
    """
    #z = [ 0.7, 1. , 1.35, 1.7, 2., 2.3, 2.6, 3., 3.5, 4.]
    #z =[0.2, 0.45, 0.7, 1.0, 1.3, 1.6, 1.85, 2.15, 2.35, 2.55, 2.85, 3.15, 3.5]
    z = [0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
    #z=[0.1,0.3,0.4,0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
    #z_mid =[0.22, 0.35, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
    z_mid =[0.32, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
    print 'this is num ', ind
    z_min  = z[ind -1]
    z_max  = z[ind]
    print z_min,z_max
    z_mean = numpy.mean((z_min,z_max))
    z_mean = z_mid[ind -1]
    return z_min,z_max, z_mean
def get_num(z_i):	
   '''
   returns the index number corrisponding to z
   '''
   z_i=round(z_i,2)
   z = [0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
   #z=[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,1.9,2.2,2.5, 2.8,3.2,3.6,4.0]
   #z_mid =list(numpy.around(medianArray(z),decimals=2))
   z_mid =[0.32, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
   #z_mid =[0.22, 0.35, 0.53, 0.7, 0.89, 1.11, 1.45, 1.75, 2.23, 2.83, 3.44]
   #print z_mid, z_i
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
Returns the comoving radial distance in Mpc
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
    sbins = fbins/(pi*4 * (dl*3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX+1))) #sbins = fbins/(pi*4 * (dl* 3.08e22)**2 * (1 + z)**((LF_SPEC_INDEX)+1))
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
        V(Vmax) - Mpc^-3
        dl      - Mpc
    Returns:
        dnds - 
    """

    L = numpy.array(lbins)
    phi_m = LF/(numpy.log(10)*L) #dLogL
    #phi_m = LF*(numpy.log(10**0.4))/L
    dndl = phi_m * V
    dnds = dndl*dlds
    #print dnds*1e-26
    return dnds*1e-26

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
def dNdS_LF(S, z, dlds, Vmax, dl, mag, params=None, paramsList=None, inta=None, area=None, family=None):
    """
    Source count model
    S is in Jy
    Cosmology is set globally
    """
    
    L=get_Lbins([S],z,dl)[0]*(1.4/3.)**(-0.7) #doing everything at 1.4GHz
    area = area*sqDeg2sr
    phi = LF(S, z, dlds,Vmax,dl,params,paramsList,inta,area,family,L=L)
    return LFtodnds(L*(1.4/3.)**(0.7),phi,dlds,Vmax,dl,mag)
   
def LF(S, z, dlds, Vmax, dl, params=None, paramsList=None, inta=None, area=None, family=None,L=None):

    #Lmin=Lmax=Lnorm=Lstar=Lslope=Lslope2=Lzevol=-99.0
    Lmin=params[paramsList.index('LMIN')]
    Lmax=params[paramsList.index('LMAX')]
    if family in ['LFsch','LFdpl_dpl','LFdpl_dpl_z','LFdpl_pl','LFpl_dpl','LFlognorm_dpl','LFpl_lognorm']:
        Lnorm=params[paramsList.index('LNORM')]
        Lstar=params[paramsList.index('LSTAR')]
        Lslope=params[paramsList.index('LSLOPE')]
        #Lzevol=params[paramsList.index('LZEVOL')]
    if family in ['LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFdpl_dpl_z','LFlognorm_dpl']:
        Lslope2=params[paramsList.index('LSLOPE2')]
    if family in ['LFlognorm','LFpl','LFpl_dpl','LFdpl','LFdpl_pl', 'LFlognorm_dpl','LFdpl_dpl','LFdpl_dpl_z','LFpl_lognorm']:
        Lnorm_2=params[paramsList.index('LNORM_2')]
        Lstar_2=params[paramsList.index('LSTAR_2')]
        Lslope_2=params[paramsList.index('LSLOPE_2')]

    if family in ['LFlognorm_dpl','LFdpl_pl','LFpl_dpl','LFdpl_dpl','LFdpl_dpl_z','LFdpl_dpl_z','LFpl_lognorm']:
        Lmin2=params[paramsList.index('LMIN2')]
        Lmax2=params[paramsList.index('LMAX2')]
        Lmin2,Lmax2 =numpy.power(10,[Lmin2,Lmax2])


    if family in ['LFdpl', 'LFdpl_dpl','LFdpl_dpl_z']:
        Lslope2_2=params[paramsList.index('LSLOPE2_2')]

    if family in ['LFdpl_dpl_z','LFevol_dpl','LFevol_logn','LFevol_logn_all','LFevol_logn_all_L','LFevol_dpl_L','LFevol_logn_L','LFevol_phi_logn']:
	alpha_agn=params[paramsList.index('A_agn')]
	alpha_SF=params[paramsList.index('A_SF')]
	beta_agn=params[paramsList.index('B_agn')]
	beta_SF=params[paramsList.index('B_SF')]

    if family in ['LFevol_dpl_s','LFevol_logn_s']:
        alpha_SF=params[paramsList.index('A_SF')]
        beta_SF=params[paramsList.index('B_SF')]

    if family in ['LFevol_dpl_a','LFevol_logn_a']:
        alpha_agn=params[paramsList.index('A_agn')]
        beta_agn=params[paramsList.index('B_agn')]

    if family in ['LFevol_logn_mat','LFevol_phi_logn_mat','LFevol_logn_el','LFevol_logn_sigma','LFevol_logn_lmin','LFevol_logn_lnorm']:
        alpha_SF=params[paramsList.index('A_SF')]
        alpha_agn=params[paramsList.index('A_agn')]

    if family in ['LFevol_logn_lmin']:
	Lmin2=params[paramsList.index('LMIN2')]
    
    if family in ['LFevol_logn_lnorm','LFevol_logn_all','LFevol_logn_all_L']:
        Lnorm=params[paramsList.index('LNORM')]
	Lstar=params[paramsList.index('LSTAR')]



    if family in ['LFevol_logn_all','LFevol_logn_slope','LFevol_logn_all_L']:
        alpha_SF=params[paramsList.index('A_SF')]
        alpha_agn=params[paramsList.index('A_agn')]
	Lslope2_2=params[paramsList.index('LSLOPE2_2')]

    if  '_phi' in family:
        alpha_D=params[paramsList.index('A_D')]
        beta_D=params[paramsList.index('B_D')]
	

    if family in ['LFevol_dpl_L','LFevol_logn_L','LFevol_logn_all_L']:
	num=get_num(z)
	#print z, num
        if num>0:
	  Lmin=params[paramsList.index('LMIN_%d'%num)]
	  #print z, num, Lmin, 'using LMIN_%d'%num
	    
    if family in ['LFlognorm','LFlognorm_dpl','LFpl_lognorm','LFevol_logn_sigma','LFevol_logn_all','LFevol_logn_all_L']:
        Lsigma = params[paramsList.index('LSIGMA')]
    #print params
    if L==None:
    	L=get_Lbins([S],z,dl)[0]*(1.4/3.)**(-.7)
	print 'come hither'
        
    #Lmin=params[paramsList.index('LMIN')]
    #Lmax=params[paramsList.index('LMAX')]
    Lmin,Lmax =numpy.power(10,[Lmin,Lmax])
    #Lmax =numpy.power(10,[Lmax])
    #Smax=get_sbins([Lmax],z,dl)
    #Lmin = get_Lbins([Lmin/1e6],z,dl)[0]*(1.4/3.)**(-.7)
    #print Lmin,L,Lmax,S#,Lnorm_2,Lstar_2,Lslope_2#,Lslope2
    #print Lmin,L,Lmax,S,Lnorm,Lstar,Lslope,Lslope2,Lnorm_2,Lstar_2,Lslope_2,Lsigma
    #print LSTAR_MIN,LSTAR_MAX
    #print family
    #sys.exit()
    #print z
    if Lmin < L < Lmax:
        if inta is not None:
            return float(inta(S))
        else:
            #print  Lmin,Lmax, L
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

                if Lmin2<L<Lmax:
		    Lr=10**(lmin2-Lstar)
		    Lr_2=10**(lmin2 - Lstar_2)
		    #Lnorm_2=numpy.log10(10**Lnorm*(Lr_2)**Lslope_2/((Lr)**Lslope + (Lr)**Lslope2))
		    #print 'Lnorm_2', 'Lr', 'Lr_2'
		    #print Lnorm_2, Lr, Lr_2
                    phi_1=powerlaw(L,Lstar_2,Lslope_2,Lnorm_2)
                if Lmin<L<Lmax2:
                    phi_2=doublepowerlaw(L,Lstar,Lslope,Lslope2,Lnorm)
                phi = phi_1+phi_2
		#print phi
		#sys.exit()


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
                        #print '-',
                            return -1.0e99
                        phi_1,phi_2 =0,0
                        #L=L*(3/1.4)**(-.7)
                        if Lmin2<L<Lmax:
                            L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
			    phi_1=doublepowerlaw(L_1,Lstar_2,Lslope_2,Lslope2_2,Lnorm_2)
                        if Lmin<L<Lmax2:
                            L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
		            phi_2 =doublepowerlaw(L_2,Lstar,Lslope,Lslope2,Lnorm)
                        phi = phi_1+phi_2
            elif family=='LFevol_dpl':
                        L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
			L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
			#print numpy.log10([L,L_1,L_2]), z
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=doublepowerlaw(L_2,22.41 ,2.4,0.43,-3.50)
                        phi = phi_1+phi_2
            elif family=='LFevol_dpl_s':
                        L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                        phi_2=doublepowerlaw(L_2,22.41 ,2.4,0.43,-3.50)
                        phi = phi_2


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
			if L_1< Lmin2:phi_1=0
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




            elif family in ['LFevol_logn', 'LFevol_logn_L'] :
                        L_1 = L/(1 + z)**(alpha_agn+z*beta_agn)
                        L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                        #print numpy.log10([L,L_1,L_2]), z
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63,numpy.log10(3.55e-3))
                        phi = phi_1+ phi_2

            elif family in ['LFevol_phi_logn'] :
                        L_1 = L/(1 + z)**(alpha_agn)
                        L_2 = L/(1 + z)**(alpha_SF)
                        phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
                        phi_2_evol = (1+z)**(alpha_D)
                        phi = phi_1+ phi_2_evol*phi_2

                        phi = phi_1+phi_2
	    elif family=='LFevol_logn_mat':
		            L_1 = L/(10**(alpha_agn))
		            L_2 = L/(10**(alpha_SF))		            
		            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
		            phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
			    #phi_2=lognormpl(L_2,21.16,1.22,0.63, -2.6)
		            phi = phi_1+phi_2
	    elif family=='LFevol_logn_el':
                            L_1 = L/(1+z)**(alpha_agn)
                            L_2 = L/(1+z)**(alpha_SF)
                            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                            #phi_2=lognormpl(L_2,20.8,1.05,0.63, -2.13) #z_6e
			    phi_2=lognormpl(L_2,21.36,1.15,0.5, -2.36) #z_6f
			    #phi_2=lognormpl(L_2,21.3,1.14,0.51, -2.32) #z_6g
			    phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
                            phi = phi_1+phi_2

            elif family=='LFevol_phi_logn_mat':
                            #L_1 = L/(10**(alpha_agn))
                            #L_2 = L/(10**(alpha_SF))
			    L_1 = L/(1+z)**(alpha_agn)
			    L_2 = L/(1+z)**(alpha_SF)
                            phi_1=doublepowerlaw(L_1,24.59,1.27,0.49,numpy.log10(2.5*10**(-5.5)))
                            phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
			    phi_2_evol = (1+z)**(alpha_D)
                            phi = phi_1+ phi_2_evol*phi_2


            elif family=='LFevol_logn_s':
                        L_2 = L/(1 + z)**(alpha_SF+z*beta_SF)
                        phi_2=lognormpl(L_2,numpy.log10(1.85e21),1.22,0.63, numpy.log10(3.55e-3))
                        phi = phi_2




            elif family=='LFlognorm':
                phi=lognormpl(L,Lstar_2,Lslope_2,Lsigma,Lnorm_2)
            elif family=='LFlognorm_dpl':
                if Lslope < Lslope2:
                    #print '-',
                    return -1.0e99
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
       inta=None,area=None,family=None,NOISE=None):
    """
    The integrand is the product of dN/dS_LF (count model) and G (the Gaussian)
    Schechter: ['LNORM','LSTAR','LSLOPE','LMIN','LMAX','LZEVOL'] + ['noise']
    Double PL: ['LNORM','LSTAR','LSLOPE','LSLOPE2','LMIN','LMAX','LZEVOL'] + ['noise']
    """
    #print paramsList
    #print params[0],params[1],params[2],params[3]
    Lmin=10**params[paramsList.index('LMIN')]
    #Lmin=params[paramsList.index('LMIN')]
    Lmax=10**params[paramsList.index('LMAX')]
    #[Smin,Smax]=get_sbins([Lmin,Lmax],redshift,dl)
    sigma=params[paramsList.index('noise')]
    if family in ['LFevol_dpl_L','LFevol_logn_L']:
        num=get_num(redshift)
        #print z, num
        if num>0:
          Lmin=10**params[paramsList.index('LMIN_%d'%num)]
	  #print 'bin', redshift, Lmin, num, 'Lmin_%d'%num
    [Smin,Smax]=get_sbins([Lmin,Lmax],redshift,dl)
    #Smax=get_sbins([Lmax],redshift,dl)
    #print redshift, Smin
    #print Smin,Smax,numpy.log10(get_Lbins([Lmin],redshift,get_dl(redshift))), numpy.log10(Lmax),
    if not floatNoise:
    	sigma = NOISE


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
                family=None,dump=None,verbose=False,inta=None,dlds = None,redshift=None,Vmax=None,dl=None,mag=None,noise=0):

    """
    For LF,
    function to calculate mock data for a given power law or interpolation object
    """

    #dl=cosmocalc(redshift,H0=Ho,WM=wm)['DL_Mpc']
    nbins=len(bins)
    II = numpy.zeros(nbins-1)
    #print paramsList
    #print params
    #print params[0], params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10],params[11],params[12],params[13],params[14],params[15]
    for ibin in xrange(nbins-1):
        sqDeg2srr=sqDeg2sr
        #print bins[ibin],bins[ibin+1]
        #sqDeg2srr=1.0
        #print area,sqDeg2sr
        IIi = abs(IL(dNdS_LF,redshift,dlds,Vmax,dl,mag,params,paramsList,bins[ibin],bins[ibin+1],\
                    inta=inta,area=sqDeg2srr*area,family=family,NOISE=noise))
        #II[ibin]=abs(IL(dNdS_LF,redshift,dlds,Vmax,dl,params,paramsList,bins[ibin],bins[ibin+1],\
        #            inta=inta,area=sqDeg2srr*area,family=family))
        II[ibin] = IIi
        #print  bins[ibin],bins[ibin+1], IIi

    #print II
    #sys.exit()
    return II

#-------------------------------------------------------------------------------
