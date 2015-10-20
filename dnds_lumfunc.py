'''
This is dnds_lumfunc.
Eliab Malefahlo
October 2015
 
This function converts dnds into LF
Usage:
Its not a stand-alone function(no main())

It takes in 
dnds, sbins(uJy), z_min,z_max 
call get_lumfunc(dnds,sbins,z_min,z_max)
suggested sbins=numpy.logspace(1,2.88,20)

returns
LF,Lbins
'''
import numpy
import os,sys,math,shutil
import importlib
from pylab import*
from matplotlib.ticker import AutoMinorLocator
from cosmocalc import cosmocalc
#import pymultinest
#from bayestackClasses import countModel
#from utils import sqDeg2sr,fetchStats
#from plat import*
#from countUtils import calculateDnByDs,medianArray

'''
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

'''
    
Ho=71
wm=0.27
fmin = 1e-29
specindx = -0.7
SURVEY_AREA = 2.62265546062 #BOSS area. avoid importing settings file
sqDeg2sr=4.0*pi*pi/129600.0

def get_Vmax(zlo,zup):
	z  = zup
	V1 = cosmocalc(z,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
	V2 = cosmocalc(zlo,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
	area = SURVEY_AREA*sqDeg2sr
	vmax = area*(V1-V2)/(4*pi)
	#print 'vmax',vmax,V1,V2,zup,zlo,SURVEY_AREA,sqDeg2sr
	return vmax
	
def get_dnds(dnds_ecl,sbins):
	dnds=[]
	for i in range(len(sbins)):
		s2_5 = 1e-26*(sbins[i]*1e-6)**2.5
		print s2_5,dnds_ecl[i],sqDeg2sr*SURVEY_AREA
		dndsi=(dnds_ecl[i]SURVEY_AREA)/(s2_5)
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
	bin_norm = numpy.zeros(len(sbins)-1)
	err_bin  = numpy.zeros(len(sbins)-1)
	rho_m 	 = numpy.zeros(len(sbins)-1)
	
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
	print o_Vmax,sqDeg2sr*SURVEY_AREA
	dnds = get_dnds(dnds_ecl,medianArray(sbins))
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



