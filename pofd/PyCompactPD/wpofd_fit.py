#!/usr/bin/env python

import os,sys
from os.path import expanduser
home = expanduser("~")

pp='.local/lib/python2.7/site-packages/'
sys.path.insert(0,os.path.join(home,pp))
os.environ['PYTHONPATH']=os.path.join(home,pp)

import numpy as np
from _wpofd6 import ffi,lib

import countUtils
from utils import sqDeg2sr,strictly_increasing

from mpi4py import MPI
import pymultinest
from priors import Priors

VERSION='_v2' # or ''

#-------------------------------------------------------------------------------

def poissonLhood(data,realisation,silent=True):
    if not silent:
        for i in range(len(data)):
            print i,data[i],realisation[i]
    kk=data[np.where(realisation > 0) and np.where(data>0)];
    iii=realisation[np.where(realisation > 0) and np.where(data>0)]
    loglike = (kk*np.log(iii) + kk - kk*np.log(kk) - iii).sum()
    return loglike

#-------------------------------------------------------------------------------

def dNByDsToArray(params,paramNames,array,gridlength=None):

    """
    Given a dn_by_ds function, populate an array
    """

    if gridlength is None:
        gridlength=len(array)

    #array=np.linspace(-20.0,100.0,22)

#    print 'aa',array # !!! array is not used for ppl
    draws,func=countUtils.simulate('ppl',params,paramNames,\
                                  array,seed=1234,N=1,noise=1.0e-10,dump=None,\
                                  area=1.0,gridlength=gridlength,\
                                  output=None,verbose=False,return_vals=True)

  #  bins=np.linspace(array[0],array[-1],gridlength)
  #  counts=np.histogram(draws,bins=bins)[0]
  #  dn_by_ds,dn_by_ds_eucl,dn_by_ds_errs,dn_by_ds_b,dn_by_ds_b_errs=\
  #    countUtils.calculateDnByDs(1.0e-6*bins,counts,idl_style=False,return_all=True)
#    print dnds_vals.size,dnds_vals
    dnds_vals=[func(1.0e6*atflux) for atflux in array]
    dn_by_ds=dnds_vals
    #print 'wwww',dnds_vals
    return dnds_vals,dn_by_ds

#-------------------------------------------------------------------------------


pri=Priors()
def logprior_pofd(cube,ndim,nparams):

    # ['C','S0','S1','S2','a0','a1']
    cube[0]=pri.GeneralPrior(cube[0],'LOG',1.0e10,1.0e20)
    cube[1]=pri.GeneralPrior(cube[1],'U',0.1,20.0)
    cube[2]=pri.GeneralPrior(cube[2],'U',0.1,100.0)
    cube[3]=pri.GeneralPrior(cube[3],'U',20.0,100.0)
    cube[4]=pri.GeneralPrior(cube[4],'U',-2.5,-0.1)
    cube[5]=pri.GeneralPrior(cube[5],'U',-2.5,-0.1)
    #print [cube[i] for i in range(nparams)]

    return


#-------------------------------------------------------------------------------

inited=False
log=open(os.path.join('chains_test6','draws.txt'),'w')
kk=0
def loglike_pofd(cube,ndim,nparams):
    """
    """

    global Nbins,dd,data,DataArray,seed,c,inited,log,kk

    if not inited:
        datadir='%s/Dropbox/bayestack/pofd/CompactPD%s/Debug/out/'%(home,VERSION)
        histof='histo.txt'
        histof=os.path.join(datadir,histof)
        data=np.loadtxt(histof)
        Nbins=data.shape[0]
        dd=np.ascontiguousarray(data.T)
        DataArray = ffi.cast("double *",dd.ctypes.data)
        #z=np.zeros(Nbins)
        #result=ffi.cast("double
        #*",np.ascontiguousarray(z).ctypes.data)
        #ParamsArray = ffi.new("struct PD_params *")
        c=[2.0e15,0.08,50.0,85.0,-1.80,-2.1]
        c=[1.39619466740548425e15,1.43870273719837827,71.0828702687349079,96.4907460253095053,\
           -2.49523428745622500,-02.19798785182257417]
        inited=True

    paramList=['C','S0','S1','S2','a0','a1']
    if not strictly_increasing([cube[i] for i in range(ndim) if paramList[i].startswith('S')]):
            #print [cube[i] for i in range(ndim)]                               
            print '+',
            return -1.0e99
    else:
        #print [cube[i] for i in range(ndim)]
        pass
    #print data[-1,:]
    #fluxes=np.array([2.0e-8,50.0e-6,85.0e-6])
    fluxes=np.linspace(2.0e-8,85.0e-6,8)
    gl=100

   # cube=range(len(paramsList))
    #print 'cc',[cube[i] for i in range(ndim)],dnds_vals

    #c[0]*=10.0
    #cube=[c[ic] for ic in range(ndim)]
    dnds_vals,dn_by_ds=dNByDsToArray([cube[i] for i in range(ndim)],\
                                     paramList,fluxes,gridlength=gl)
    #print 'vv',dnds_vals
    x=np.log10(np.array(fluxes))
    power=0.0
    y=np.log10(np.array(dnds_vals)*np.power(np.power(10,x),power))
    #print 'xy',y
    model=np.array([x,y])
    model[-1,0]=model[-1,1]
    #print 'wwxx',model[-1,:]
    #np.savetxt('xy.txt',model.T)
    #(model,ParamsArray)=pars
    strr=str(kk)+' cc '+str([cube[i] for i in range(ndim)])+' '+np.array_str(model)
    print strr
    log.write(strr+'\n')
    log.flush()
    #xy=model
    #print 'zz',model[-1,:]#model.min(),model.max()

 #   i_p=ffi.cast("double *",np.ascontiguousarray(model).ctypes.data)
    i_l=model.shape[-1]
    #print 'il',i_l


    #print type(model),model.shape


    #i_l=10
    # Model is changing
    # i_p is not...
    # You can change model and i_p doesn't change
    #model[-1,:]=2.0*np.ones((1,10))
   # model[-1,:]=np.random.uniform(8.314,16.1737,i_l)
    i_p=ffi.cast("double *",np.ascontiguousarray(model).ctypes.data)
    #print 'nn0',model[-1,:]
    ParamsDict={\
        'd_max':1.0e-3,
        'd_min':10**-7.32,
        'source_max':8.5e-5,
        'source_min':10**-7.32,
        'PSFresultionFWHM':6.0,
        'pixelsize':1.0,
        'sigma_noise':17.0e-6,
        'interplot_length':i_l,
        'interplot_pointer':i_p}#ffi.cast("double *",np.ascontiguousarray(model).ctypes.data)}
    ParamsArray = ffi.new("struct PD_params *",ParamsDict)

    m=model
#    #np.random.seed(seed=int(np.random.uniform(0,100,1)))
#    #m[-1,:]=np.linspace(1.0,10.0,10)
#    #m[-1,:]=np.random.uniform(0.0,1.0,10)
#    #print m

#    ParamsArray.interplot_pointer=ffi.cast("double *",\
#                                           np.ascontiguousarray(m).ctypes.data)
#    print ParamsArray.interplot_pointer[0]
#    p=[ParamsArray.interplot_pointer[i+i_l] for i in range(i_l)]
#    print 'pp',p
#    sys.exit(0)

    z=np.zeros(Nbins)
    result=ffi.cast("double *",np.ascontiguousarray(z).ctypes.data)
#    print 'ii',i_p[0]
    #print 'pp',ParamsDict['interplot_pointer'][101]
    r=lib.CompactPD_LH(Nbins,DataArray,result,ParamsArray)
    #print 'rrr',result[101]
    #print 'rrxx',r
    realisation=np.array([result[i] for i in range(Nbins)])
    #print 'yy1',z.sum(),r
    #print 'ffxx',realisation
    #np.savetxt('realn.txt',realisation)
    #sys.exit(0)
    loglike=poissonLhood(data[:,-1],np.power(10,realisation),silent=True)

    lk=str(kk)+' ll1 '+str(loglike)
    print lk
    log.write(lk+'\n')
    log.flush()
    kk+=1
#    print kk,'ll1',loglike
    #print m[0,:]
    #print 'junk',m
    #m=np.array([[-7.32,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9],\
    #            [16.1737,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,8.314]])
 #   m[-1,:]=np.random.uniform(8.314,16.1737,i_l)
    #print 'nn1',m[-1,:]
 #   m2=m
#    m2[-1,:]=m[-1,:]*2.0
 #   m2[-1,:]=np.random.uniform(8.314,16.1737,i_l)
    #print 'nn2',m2[-1,:]
#    ParamsArray.interplot_pointer=ffi.cast("double *",\
#                                           np.ascontiguousarray(m2).ctypes.data)
#    p=[ParamsArray.interplot_pointer[i+i_l] for i in range(i_l)]
#    print 'pp',p
    del ParamsArray,ParamsDict,realisation,result,i_l,i_p
    #sys.exit(0)
   # z=np.ones(Nbins)
    # Even though ParamsArray.interplot_pointer is changing, result(2) is not
   # result2=ffi.cast("double *",np.ascontiguousarray(z).ctypes.data)
   # r=lib.CompactPD_LH(Nbins,DataArray,result2,ParamsArray)

   # z=np.array([result2[i] for i in range(Nbins)])
   # print 'yy2',z.sum(),r
    #sys.exit(0)
   # realisation=z
    #del result,ParamsArray,i_p
   # print 'ff',realisation

   # loglike2=poissonLhood(data[:,-1],np.power(10,realisation),silent=False)
   # print 'll2',loglike2
    #sys.exit(0)
    return loglike

#-------------------------------------------------------------------------------

def main():

    global nparams,ndim
    
    #Npix=data[:,-1].sum()
    #Nbins=data.shape[0]
    #dd=np.ascontiguousarray(data.T)
    #DataArray = ffi.cast("double *",dd.ctypes.data)

    if '2' in VERSION:
        """
        Takes two arrays x, y for the model coordinates
        """
        #x=[-7.32,-6.7,-6.3,-5.533,-4.766,-4.0,-3.25,-1.9]
        #y=[16.1737,15.0620,14.4342,13.4563,12.1586,10.2764,8.55,6.314]

        # These are very roughly the SKADS counts
        #fluxes=[2.0e-8,17.0e-6,85.0e-6]
        #dnds_eucl_per_sr=[5.0e-3,1.0,3.0]

        paramList=['C','S0','S1','S2','a0','a1']
        params=[4000.0,0.08,50.0,85.0,-1.80,-2.1]
        params=[1780.0,25.0,49.0,74.5,-0.129,-0.129]

        fluxes=np.array([2.0e-8,50.0e-6,85.0e-6])
        fluxes=np.linspace(2.0e-8,85.0e-6,10)
        print fluxes
        gl=10
        dnds_vals,dn_by_ds=dNByDsToArray(params,paramList,1.0e6*fluxes,gridlength=gl)
        fluxes=np.linspace(fluxes[0],fluxes[-1],gl)
        #print dnds_vals,dn_by_ds*sqDeg2sr
        #sys.exit(0)

        # Some housekeeping
        x=np.log10(np.array(fluxes))
        power=0.0
        y=np.log10(np.array(dnds_vals)*np.power(np.power(10,x),power))
        #y=dnds_vals
        #print dnds_vals
        #y=dn_by_ds
        #print x,x.size
        #print y,y.size
        print y
        print x.shape,y.shape

        model=np.array([x,y])
        np.savetxt('xy.txt',model.T)
        #sys.exit(0)
        i_p=ffi.cast("double *",np.ascontiguousarray(model).ctypes.data)
        i_l=model.shape[-1]

#        ParamsArray = ffi.new("struct PD_params *",{\
#        'd_max':1.0e-3,
#        'd_min':10**-7.32,
#        'source_max':8.5e-5,
#        'source_min':10**-7.32,
#        'PSFresultionFWHM':6.0,
#        'pixelsize':1.0,
#        'sigma_noise':17.0e-6,
#        'interplot_length':i_l,
#        'interplot_pointer':i_p})

        #z=np.zeros(Nbins)
        #result=ffi.cast("double *",np.ascontiguousarray(z).ctypes.data)
        #r=lib.CompactPD_LH(Nbins,DataArray,result,ParamsArray)
        #realisation=np.array([result[i] for i in range(Nbins)])
        #loglike=poissonLhood(data[:,-1],np.power(10,realisation),silent=True)
        #print 'loglike = %f [unchecked]\n'%loglike
#        pars=(model,ParamsArray)
#        print '# ibin data realn'
        #loglike=loglike_pofd(pars)
        #print 'loglike = %f [unchecked]\n'%loglike

        nparams=ndim=len(params)
        pymultinest.run(loglike_pofd,logprior_pofd,nparams,verbose=True,multimodal=False,\
                        write_output=True,n_live_points=500,resume=True,mode_tolerance=-1e90,
                        seed=1234,max_iter=0,importance_nested_sampling=False,\
                        outputfiles_basename='chains_test6/1-',init_MPI=False)

    else:
        """
        Takes a data array but uses fixed model (for testing)
        """
        ParamsArray = ffi.new("struct PD_params *",{\
        'd_max':1.0e-3,
        'd_min':np.power(10,-7.32),
        'source_max':8.5e-5,
        'source_min':np.power(10,-7.32),
        'last_interplot_log10x':-7.32,
        'last_interplot_log10y':16.1737,
        'PSFresultionFWHM':6.0,
        'pixelsize':1.0,
        'sigma_noise':17.0e-6})

        result=lib.CompactPD_LH(Nbins,DataArray,ParamsArray)
        print '\nresult: %f [OK]\n' % result
        chisq_true=1553478.027942; tol=1.0e-3
        assert(result-chisq_true<tol), '***FAILURE: chisq is incorrect!!!'

    return 0

#-------------------------------------------------------------------------------

def main2():

    import numpy as np
    from cffi import FFI

    ffi = FFI()
    ffi.cdef("void copy(float *in, float *out, int len);")
    C = ffi.dlopen("ccode.dll")

    a = 42*np.ones(16, dtype=np.float32)
    b = np.zeros_like(a)
    pa = ffi.cast("float *", a.ctypes.data)
    pb = ffi.cast("float *", b.ctypes.data)

    C.copy(pa, pb, len(a))
    print b

    return 0

#-------------------------------------------------------------------------------

if __name__=='__main__':

    ret=main()
    sys.exit(ret)
