'''
This is dnds_lumfunc.
Eliab Malefahlo
October 2015
 
This function converts dnds into LF
Usage:

Arguments
dnds  : The reconstructed source counts (from reconstruct)
sbins : The flux(uJy) bins that coresponds to dnds
z_min : minimum redshift of slice
z_max : max redshift of slice

returns 
rho_m : the normalized luminosity function (Log[Mpc^-3 mag^-1])
Lbins : Log Luminosity bins corresponding to rho_m (Log[L])

dnds = [4.78118760e-05,   6.33563041e-05,   8.00873144e-05,   8.50948383e-05]
sb = [ 45.,     55. ,    65.,     75.]
dnds_lumfunc.get_lumfunc(dnds,sb,1.8,2.5)

output
(array([ -7.60008913e+00,  -7.52893349e+00,  -7.46877069e+00,
         1.21271141e-33]),
 array([ 24.35010723,  24.43725741,  24.50980808,  24.57195598]))

'''
import numpy
import os,sys,math,shutil
import importlib
from pylab import*
from matplotlib.ticker import AutoMinorLocator
from cosmocalc import cosmocalc
import pymultinest
from bayestackClasses import countModel
from utils import sqDeg2sr,fetchStats
from plat import*
from countUtils import calculateDnByDs,medianArray


param_file=sys.argv[-1]
setf='%s.bayestack_settings' % param_file

print '%s is using %s' % (__name__,setf)

set_module=importlib.import_module(setf)
globals().update(set_module.__dict__)

expt=countModel(modelFamily,nlaws,setf,[dataset],floatNoise)
#print expt.data
print len(expt.bins)
print expt.bins
# Fetch best-fit parameters and calculate best-fit line
plotTruth=dict((name,-99.0) for name in expt.parameters)
stats=fetchStats(outdir,expt.parameters,plotTruth)
SMIN_MAP=stats['S0'][0]
SMIN_MAP_UPPER= stats['S0'][-2]
SMIN_MAP_LOWER=stats['S0'][-1]

s=loadtxt('%s/recon_stats.txt'%outdir)

z = [ 0.7, 1 , 1.35, 1.7, 2, 2.3, 2.6, 3, 3.5, 4]
if 'sdss' in datafile:
	ind = int(datafile[-5])
	z_min= z[ind -1]
	z_max= z[ind]
	print z_min,z_max
    

def get_Vmax(zlo,zup):
	z  = zup
	V1 = cosmocalc(z,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
	V2 = cosmocalc(zlo,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
	area = SURVEY_AREA*sqDeg2sr
	vmax = area*(V1-V2)/(4*pi)
	#print 'vmax',vmax,V1,V2,zup,zlo,SURVEY_AREA*sqDeg2sr,sqDeg2sr
	return vmax
	
def get_dnds(dnds_ecl,sbins):
	dnds=[]
	for i in range(len(sbins)):
		s2_5 = 1e-26*(sbins[i]*1e-6)**2.5
		print s2_5,dnds_ecl[i],SURVEY_AREA*sqDeg2sr
		dndsi=(dnds_ecl[i]*SURVEY_AREA*sqDeg2sr)/(s2_5)
		dnds.append(dndsi)
	return numpy.array(dnds)

def get_dndl(dnds,z,dl):
	dsdl = 1./(math.pi*4 * (dl* 3.08e22)**2 * (1 + z)**((-0.7)+1))
	print z
	print 'dsdl'
	print dsdl
	return dnds*dsdl
	
def get_Lbins(sbins,z,dl):
	Lbins = [math.pi*4 *(s*1e-32)* (dl* 3.08e22)**2 * (1 + z)**((-0.7)+1) for s in sbins]
	return Lbins
	
def get_lumfunc(dnds_ecl,sbins,z_min,z_max):
	rho_m 	 = numpy.zeros(len(sbins)-1)
	dnds_ecl = numpy.array(dnds_ecl)
	sbins    = numpy.array(sbins)
	print numpy.where(sbins>0)
	
	dnds_ecl = dnds_ecl[numpy.where(sbins>0)]
	sbins    = sbins[numpy.where(sbins>0)]
	print 'recon'
	print dnds_ecl
	
	#print len(dnds_ecl), len(sbins)
	#print sbins
	z  = mean((z_min,z_max))
	dl = cosmocalc(z,H0=Ho,WM = wm)['DL_Mpc']
	
	Vmax=get_Vmax(z_min,z_max)
	o_Vmax    =  1./(Vmax)
	print '1/vmax'
	print o_Vmax,SURVEY_AREA*sqDeg2sr
	dnds = get_dnds(dnds_ecl,sbins)
	print 'dsdn'
	print dnds
	dndl = get_dndl(dnds, z,dl)
	rho_m 	= dndl*o_Vmax
	print 'dndl'
	print dndl
	print 'rho'
	print rho_m
	#err_bin = dndl*(o_Vmax**2)
	L = get_Lbins(sbins,z,dl) 
	Lbins =log10(L)
	for i in range(len(sbins)-1):
		dm         = (Lbins[i+1] - Lbins[i])/0.4
		rho		   = rho_m[i]/(dm)*L[i]
		rho_m[i]   = log10(rho)
		#lin_rho[i] = rho

	#err_bin = err_bin**0.5
	#error = [abs(log10(rho + err)) - abs(log10(rho)) for(err,rho) in zip(err_bin,lin_rho)]
	
	print rho_m,Lbins
	return rho_m,Lbins



