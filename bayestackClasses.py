"""
Support classes for bayestack etc.

Jonathan Zwart
May 2015

Especially countModel, e.g.

expt=countModel(modelFamily,nlaws,settingsf,dataset,floatNoise)

then exposes

expt.parameters
expt.data
expt.bins
expt.nbins
expt.binsMedian
expt.realise()
expt.logprior()
expt.loglike()
expt.dn_by_ds()
expt.confusionNoiseSquared()

etc.

"""


import os
import importlib
import numpy
import itertools
from scipy import integrate
from scipy.special import erf
from priors import Priors
import countUtils
import polnUtils
import lumfuncUtils
from utils import sqDeg2sr,beamFac,sqrtTwo,strictly_increasing,poissonLhood,medianArray,poissonLhoodMulti2
from cosmocalc import cosmocalc
#import cosmolopy
from collections import OrderedDict

#-------------------------------------------------------------------------------

class surveySetup(object):

    """
    Set up the survey parameters (assumed fixed)
    survey=surveySetup(whichSurvey)
    
    """

    def __init__(self,whichSurvey,datafiles,areas,noises):
        self.whichSurvey=whichSurvey
        # Start handling multiple datafiles
        self.datafiles=datafiles
        self.SURVEY_AREAS=areas
        self.SURVEY_NOISES=noises
        if len(self.datafiles)>1:
            self.multi=True
            self.fractions=[a/sum(self.SURVEY_AREAS) for a in self.SURVEY_AREAS]
        else:
            self.multi=False
            self.datafile=datafiles[0]
        self.SURVEY_AREA=areas[0]
        self.SURVEY_NOISE=noises[0]
        self.fractions=[1.0]

        if whichSurvey in ['video']:# or 'sim' in whichSurvey:
            self.HALO_MASK=11436315.0/(19354.0*19354.0)
            self.SURVEY_AREA=areas[0]*(1.0-self.HALO_MASK)# sq.deg. [Boris -> 0.97 sq. deg.]
            #self.radioSynthBeamFWHM=4.0 # pixels/upsamplingFactor
            #self.radioSynthOmegaSr=sqDeg2sr*beamFac*(self.radioSynthBeamFWHM/3600.0)**2

#-------------------------------------------------------------------------------

class dataSetup(object):
    """
    Structure for organizing multiple data sets
    Run as e.g.
    dataObject=dataSetup('sdss','sdss_manifest.txt',redshiftSlices=True)
    """

    def __init__(self,whichSurvey,manifestf,noiseZones=False,redshiftSlices=False):
        self.whichSurvey=whichSurvey
        self.noiseZones=noiseZones
        self.redshiftSlices=redshiftSlices
        self.manifestf=manifestf
        self.manifest=numpy.genfromtxt(self.manifestf)
        # Extract selected z slices
        self.manifest=self.manifest[[s-1 for s in whichRedshiftSlices],:]
        self.whichRedshiftSlices=whichRedshiftSlices
        if self.redshiftSlices:
            self.datafiles=['%s/sdss_dr12s%i.txt'%(whichSurvey,r) for r in self.whichRedshiftSlices]
            self.datafile=self.datafiles[0]
            self.numRedshiftSlices=len(whichRedshiftSlices) #self.manifest.shape[0]
            self.redshifts=self.manifest[:,3]
            self.sliceAreasTemp=self.manifest[:,4]  # Needs to be rejigged
            self.sliceNoisesTemp=self.manifest[:,5] # Needs to be rejigged
            self.zsliceData={}; self.dls={}; self.sliceAreas={}; self.sliceNoises={}
            # Rejig z information
            for r in range(self.numRedshiftSlices):
                z=self.redshifts[r]
                dataf=self.datafiles[r]
                self.zsliceData[z]=self.loadData(dataf)
                self.dls[z]=cosmocalc(z,H0=Ho,WM=wm)['DL_Mpc']
                self.sliceAreas[z]=self.sliceAreasTemp[r]
                self.sliceNoises[z]=self.sliceNoisesTemp[r]
            print 'Loaded',self.datafiles
        elif self.noiseZones:
            pass

        return

    def loadData(self,datafile):
        dataMatrix=numpy.genfromtxt(datafile)
        corrected_data=dataMatrix[:,3]*dataMatrix[:,8]
        binsDogleg=numpy.concatenate((dataMatrix[:,0],[dataMatrix[-1,1]]))
        return corrected_data,binsDogleg


#-------------------------------------------------------------------------------

class model(object):

    """
    Base class for countModel class
    A model is specified by a parameter vector
    and a function relating them
    Should be vectorized/able
    usage:
        model=model(family)
        model.eval(vals,params)

    e.g.

    import numpy
    from bayestackClasses import *

    p={'amp':2.0,'breaks':[1.0,2.0,3.5],'slopes':[-2.0,-1.5],'coeffs':[1.0,2.0,3.0]}
    v=numpy.linspace(0.0,10.0,101)

    m=model('ppl')
    m.eval(v,p)

    n=model('poly')
    n.eval(v,p)
    """

    def __init__(self,family):
        self.family=family

    def eval(self,vals,params):

        if self.family=='ppl':
            #self.func=params['amp']*cosmolopy.utils.PiecewisePowerlaw(\
            #            params['breaks'],params['slopes'],\
            #            externalval=0.0,coefficients=None)
            self.func=None # temporary
        elif self.family=='poly':
            self.func=numpy.poly1d(list(reversed(params))) # was params['coeffs']

        return self.func(vals)


#-------------------------------------------------------------------------------

class countModel(object):

    """
    usage:
    expt=countModel(kind,order,settingsf,whichSurvey,floatNoise)
    
    """

    def __init__(self,kind,order,settingsf,whichSurvey,floatNoise,doPoln=False,\
                 doRayleigh=False,doRedshiftSlices=False,mybins=None):
        # Import settings
        print 'Settings file is %s' % settingsf
        set_module=importlib.import_module(settingsf)
        globals().update(set_module.__dict__)

        # Set up model
        self.kind=kind
        self.name=kind
        self.order=order
        self.nlaws=order
        self.model=model(self.kind)
        self.floatNoise=floatNoise
        self.doPoln=doPoln
        self.doRayleigh=doRayleigh
        self.doRedshiftSlices=doRedshiftSlices
        if self.doRedshiftSlices:
            self.zDataObject=dataSetup('sdss',zmanifestf,redshiftSlices=self.doRedshiftSlices)
            self.redshifts=redshifts

        # Set up parameters for this model
        # --> This defines the order of the parameters:
        self.paramsAvail=\
                    OrderedDict([('extra',['noise']),\
                     ('breaks',['S%i'%ic for ic in xrange(self.nlaws+1)]),\
                     ('amp',['C']),\
                     ('coeffs',['p%i'%ic for ic in xrange(self.nlaws)]),\
                     ('limits',['S%i'%ic for ic in xrange(2)]),\
                     ('poles',['b%i'%ic for ic in xrange(self.nlaws)]),\

                     ('slopes',['a%i'%ic for ic in xrange(self.nlaws)]),\
                     ('schechter',['LMIN','LMAX','LNORM','LSTAR','LSLOPE','LZEVOL']),\
                     ('doublepl',['LMIN','LMAX','LNORM','LSTAR','LSLOPE','LSLOPE2','LZEVOL']),\
                     ('schechterM',['MMIN','MMAX','MNORM','MSTAR','MSLOPE','MZEVOL']),\
                     ('doubleplM',['MMIN','MMAX','MNORM','MSTAR','MSLOPE','MSLOPE2','MZEVOL'])])
        familyMap={'ppl':['breaks','slopes','amp','extra'],\
                   'poly':['limits','coeffs','extra'],\
                   'bins':['poles','extra'],\
                   'LFsch':['schechter','extra'],\
                   'LFdpl':['doublepl','extra'],\
                   'HIsch':['schechterM','extra'],\
                   'HIdpl':['doubleplM','extra']}
        self.paramsStruct=\
          [self.paramsAvail[p] for p in self.paramsAvail if p in familyMap[kind]]
        self.parameters=list(itertools.chain(*self.paramsStruct))
        self.nparams=len(self.parameters)
        print self.nparams,self.parameters
        self.currentPhysParams=-99.0*numpy.ones(self.nparams)

        # Set up priors
        self.priors=Priors()
        self.priorsDict=self.parsePriors(self.parameters,self.floatNoise)

        # Load the data and derive the bins
        self.survey=surveySetup(whichSurvey,[datafile],[SURVEY_AREA],[SURVEY_NOISE])
        self.data,self.bins=self.loadData(self.survey.datafile)
        #print self.data,self.survey.datafile
        self.nbins=len(self.bins)-1
        self.binsMedian=medianArray(self.bins)
        self.nsrc=int(self.data.sum())
        # And load any multiple data sets
        self.fdata={}; self.fbins={}; self.fnbins={}; self.fbinsMedian={}
        for df in self.survey.datafiles:
            self.fdata[df],self.fbins[df]=self.loadData(df)
            self.fnbins[df]=len(self.fbins[df])-1
            self.fbinsMedian[df]=medianArray(self.fbins[df])
        if mybins is not None:
            self.bins=mybins
            self.binsMedian=medianArray(self.bins)
            self.nbins=len(self.bins)-1

        return

    def parsePriors(self,parameters,floatNoise):
        priorsDict={}
        try: # temporary hack
            iSmax=int([i for i in parameters if i.startswith('S')][-1][-1]) # Smax
        except IndexError:
            iSmax=-1
        for p in parameters:
            if self.kind=='ppl':
                if p.startswith('C'): priorsDict[p]=[C_PRIOR,C_MIN,C_MAX] # amplitude
                elif p.startswith('S') and setBreaks: #fixed breaks
                    if nlaws > 1: priorsDict[p]=['U',S1_MIN,S1_MAX] # This only catchs S1 
               	    if nlaws > 2 and p =='S2': priorsDict[p]=['U',S2_MIN,S2_MAX]
                    if nlaws > 3 and p =='S3': priorsDict[p]=['U',S3_MIN,S3_MAX]
                elif p.startswith('S'): priorsDict[p]=['U',SMIN_MIN,SMAX_MAX] # breaks
                elif p.startswith('a'): priorsDict[p]=['U',SLOPE_MIN,SLOPE_MAX] # slopes
            elif self.kind=='poly':
                if p.startswith('p'): priorsDict[p]=[POLYCOEFF_PRIOR,POLYCOEFF_MIN,POLYCOEFF_MAX] # #coeffs
                #if p=='p0': priorsDict[p]=['DELTA',1.0,1.0]
            elif self.kind=='bins':
                if p.startswith('b'): priorsDict[p]=[POLEAMPS_PRIOR,POLEAMPS_MIN,POLEAMPS_MAX] # bins/poles/nodes

            if p.startswith('n'): # noise
                if floatNoise:
                    priorsDict[p]=['U',NOISE_MIN,NOISE_MAX]
                else:
                    priorsDict[p]=['DELTA',SURVEY_NOISE,SURVEY_NOISE]
            elif p=='S0': priorsDict[p]=['U',SMIN_MIN,SMIN_MAX] # Smin
            elif p=='S%i'%iSmax: priorsDict[p]=['U',SMAX_MIN,SMAX_MAX] # Smax

            # LFs: ['LNORM','LSTAR','LSLOPE','LMIN','LMAX','LZEVOL'] + ['LSLOPE2']
            if p=='LMIN': priorsDict['LMIN']=[LMIN_PRIOR,LMIN_MIN,LMIN_MAX]
            elif p=='LMAX': priorsDict['LMAX']=[LMAX_PRIOR,LMAX_MIN,LMAX_MAX]
            elif p=='LNORM': priorsDict['LNORM']=[LNORM_PRIOR,LNORM_MIN,LNORM_MAX]
            elif p=='LSTAR': priorsDict['LSTAR']=[LSTAR_PRIOR,LSTAR_MIN,LSTAR_MAX]
            elif p=='LSLOPE': priorsDict['LSLOPE']=[LSLOPE_PRIOR,LSLOPE_MIN,LSLOPE_MAX]
            elif p=='LSLOPE2': priorsDict['LSLOPE2']=[LSLOPE2_PRIOR,LSLOPE2_MIN,LSLOPE2_MAX]
            elif p=='LZEVOL': priorsDict['LZEVOL']=[LZEVOL_PRIOR,LZEVOL_MIN,LZEVOL_MAX]

            # MFs: ['MNORM','MSTAR','MSLOPE','MMIN','MMAX','MZEVOL'] + ['MSLOPE2']
            if p=='MMIN': priorsDict['MMIN']=[MMIN_PRIOR,MMIN_MIN,MMIN_MAX]
            elif p=='MMAX': priorsDict['MMAX']=[MMAX_PRIOR,MMAX_MIN,MMAX_MAX]
            elif p=='MNORM': priorsDict['MNORM']=[MNORM_PRIOR,MNORM_MIN,MNORM_MAX]
            elif p=='MSTAR': priorsDict['MSTAR']=[MSTAR_PRIOR,MSTAR_MIN,MSTAR_MAX]
            elif p=='MSLOPE': priorsDict['MSLOPE']=[MSLOPE_PRIOR,MSLOPE_MIN,MSLOPE_MAX]
            elif p=='MSLOPE2': priorsDict['MSLOPE2']=[MSLOPE2_PRIOR,MSLOPE2_MIN,MSLOPE2_MAX]
            elif p=='MZEVOL': priorsDict['MZEVOL']=[MZEVOL_PRIOR,MZEVOL_MIN,MZEVOL_MAX]

        return priorsDict

    def setParams(self,params):
        self.lastPhysParams=self.currentPhysParams
        self.currentPhysParams=params

    def lastParams(self):
        return self.lastPhysParams

    def loadData(self,datafile):
        dataMatrix=numpy.genfromtxt(datafile)
        corrected_data=dataMatrix[:,3]*dataMatrix[:,8]
        binsDogleg=numpy.concatenate((dataMatrix[:,0],[dataMatrix[-1,1]]))
        return corrected_data,binsDogleg

    def evaluate(self,params):
        """
        Not yet geared for LFs OR MFs
        """
#        if mybins is not None:
#            self.binsMedian=mybins
        if self.kind=='ppl':
            C=alpha=Smin=Smax=beta=S0=gamma=S1=delta=S2=-99.0
            paramsList=self.parameters
            nlaws=int(0.5*len(paramsList)-1) # infer nlaws from length of param vector
            C=params[paramsList.index('C')]
            Smin=params[paramsList.index('S0')]
            alpha=params[paramsList.index('a0')]
            if nlaws > 1: # nlaws=2,3,4
                beta=params[paramsList.index('a1')]
                S0=params[paramsList.index('S1')]
            if nlaws > 2: # nlaws=3,4
                gamma=params[paramsList.index('a2')]
                S1=params[paramsList.index('S2')]
            if nlaws > 3: # nlaws=4
                delta=params[paramsList.index('a3')]
                S2=params[paramsList.index('S3')]
            iSmax=int([i for i in paramsList if i.startswith('S')][-1][-1])
            #print 'iS',iSmax
            Smax=params[paramsList.index('S%i'%iSmax)]
            #print paramsList,iSmax,nlaws#nlaws,S0,S1,S2,Smax,C
#            print nlaws,C,alpha,beta,Smin,Smax,S0,gamma,S1,delta,S2
            #print 'C',C
#            evaluations=[countUtils.powerLawFuncWrap(nlaws,S,C,alpha,-99.0,beta,\
#                        Smin/1.0e6,Smax/1.0e6,S0/1.0e6,gamma,S1/1.0e6,delta,S2/1.0e6,\
#                        1.0) for S in self.binsMedian/1.0e6]
            power=0.0
            # C is already in sr units
            evaluations=[countUtils.powerLawFuncWrap(nlaws,S,C,alpha+power,-99.0,beta+power,\
                        -1.0e10,1.0e10,S0/1.0e6,gamma+power,S1/1.0e6,delta+power,S2/1.0e6,\
                        1.0) for S in self.binsMedian/1.0e6]

                        #self.survey.SURVEY_AREA) for S in self.binsMedian/1.0e6]
            #print *sqDeg2sr
            #sys.exit(0)
            #power=0.0
            #if nlaws==2:
            #    evaluations=[countUtils.powerLawFuncWrap(2,S,C,\
            #        alpha+power,-99.0,beta+power,\
            #        -1.0e10,1.0e10,S0/1.0e6,-99.0,-99.0,\
            #        -99.0,-99.0,1.0) for S in self.binsMedian/1.0e6]
            #elif nlaws==3:
            #    evaluations=[countUtils.powerLawFuncWrap(3,S,C,\
            #        alpha+power,-99.0,beta+power,\
            #        -1.0e10,1.0e10,S0/1.0e6,gamma+power,\
            #        S1/1.0e6,-99.0,-99.0,1.0) for S in self.binsMedian/1.0e6]
            #            #self.survey.SURVEY_AREA*sqDeg2sr) for S in self.binsMedian/1.0e6]
            #S=0.1e-6
            #print 'CC',C,alpha,beta,S0,S1,C*S**alpha
            #sys.exit(0)
            #alpha=-2.5
            #evaluations=[C*S**(alpha) for S in self.binsMedian/1.0e6]

            #print alpha,beta,gamma
            #if nlaws>2:
            #    temp=beta; beta=gamma; gamma=temp
            #print alpha
            #if nlaws==1:
            #    evaluations=[countUtils.powerLawFuncS(S,\
            #                    C,alpha,-1.0e10,1.0e10,1.0) for S in self.binsMedian/1.0e6]
            #evaluations=[100.0*countUtils.powerLawFuncS(S,C,alpha,Smin/1.0e6,Smax/1.0e6,1.0)\
            #             for S in self.binsMedian/1.0e6]
            if doPoln:
                s=1.0#*sqDeg2sr*self.survey.SURVEY_AREA#sqDeg2sr
                # ??
                #print self.survey.SURVEY_AREA
                #sys.exit(0)
                #evaluations=[e*s for e in evaluations]
                evaluations=[e*self.survey.SURVEY_AREA for e in evaluations]
                #evaluations=[e*s*sqDeg2sr/self.survey.SURVEY_AREA for e in evaluations]
            #print evaluations
            #sys.exit(0)
        elif self.kind=='poly':
            paramsList=self.parameters
            Smin=params[paramsList.index('S0')]
            Smax=params[paramsList.index('S1')]
            coeffs=[params[paramsList.index(p)] for p in paramsList if p.startswith('p')]
            S_1=1.0 # ref flux
            evaluations=[1.0 * countUtils.polyFunc(S,S_1,Smin/1.0e6,Smax/1.0e6,\
                                             coeffs) for S in self.binsMedian/1.0e6]
        elif self.kind=='LFsch':
            # XXXX
            pass
        elif self.kind=='LFdpl':
            # XXXX
            pass
        elif self.kind=='HIsch':
            # XXXX
            pass
        elif self.kind=='HIdpl':
            # XXXX
            pass

        else:
            print '***%s unsupported right now!' % self.kind
            return
        #print evaluations
        return evaluations

    def secondMoment(self,Slower,Supper):
        return self.momentN(Slower,Supper,2)

    def momentN(self,Slower,Supper,N):
        mN=integrate.quad(lambda S:S**N*self.evaluate(S),Slower,Supper)[0]
        assert(mN>=0), 'mN = %f !>=0' % mN
        return mN

    def confusionNoiseSquared(self,Slower,Supper):
        return self.survey.radioSynthOmegaSr*self.secondMoment(Slower,Supper)

    #def checkSminPriorOK(self):
    #    (C,alpha,Smin,Smax)=self.currentPhysParams
    #    self.SminPriorOK=lumfunc.checkSminPriorOK(C,alpha,Smin,Smax,\
    #                                              self.survey.SURVEY_NOISE,\
    #                                      self.survey.radioSynthOmegaSr,numerical=False)
    #    return self.SminPriorOK

    def dn_by_ds(self,return_all=None,data=None):
        if data is None: data=self.dataRealisation
        return countUtils.calculateDnByDs(self.bins/1.0e6,data,\
                                return_all=return_all)

    def erfs(self,S,Sbin,sigma):
        """vectorized"""
        Sbinlow,Sbinhigh=Sbin
        e=0.5*(erf((S-Sbinlow)/(sqrtTwo*sigma)) - erf((S-Sbinhigh)/(sqrtTwo*sigma)))
        return e

    def convertPosterior(self,draw,power):
        """
        As it stands, this is just a pass-through
        """
        for p in self.parameters:
            #if p.startswith('S') or p.startswith('n'): # uJy -> Jy
            #    self.currentPhysParams[self.parameters.index(p)] \
            #      = draw[self.parameters.index(p)]/1.0e6
            if p.startswith('a'): # x S^{power}
                self.currentPhysParams[self.parameters.index(p)] \
                  = draw[self.parameters.index(p)]
            else: # pass
                self.currentPhysParams[self.parameters.index(p)] \
                  = draw[self.parameters.index(p)]
        return self.currentPhysParams

    def realise(self,cube):
        if self.doPoln:
            #print self.survey.SURVEY_AREA
            self.dataRealisation=polnUtils.calculateP3(cube,self.parameters,\
                    family=self.kind,bins=self.bins,\
                    area=self.survey.SURVEY_AREA,doRayleigh=self.doRayleigh)
        elif self.kind in ['LFsch','LFdpl']:
            self.zRealiseObject={}
            for r in range(len(sorted(self.zDataObject.zsliceData.keys()))):
                z=self.redshifts[r]
                zbins=self.zDataObject.zsliceData[z][1]
                dl=self.zDataObject.dls[z]
                self.zRealiseObject[z]=lumfuncUtils.calculateL3(cube,self.parameters,\
                    bins=zbins,area=self.survey.SURVEY_AREA,family=self.kind,\
                    redshift=z,dl=dl)
            return self.zRealiseObject
            #self.ziter=self.redshifts[0] # Need to select z somehow!!
            #self.dataRealisation=lumfuncUtils.calculateL3(cube,self.parameters,\
            #        bins=self.bins,area=self.survey.SURVEY_AREA,family=self.kind,\
            #        redshift=self.ziter)
        elif self.kind in ['HIsch','HIdpl']:
            self.zRealiseObject={}
            for r in range(len(sorted(self.zDataObject.zsliceData.keys()))):
                z=self.redshifts[r]
                zbins=self.zDataObject.zsliceData[z][1]
                dl=self.zDataObject.dls[z]
                self.zRealiseObject[z]=hiUtils.calculateH3(cube,self.parameters,\
                    bins=zbins,area=self.survey.SURVEY_AREA,family=self.kind,\
                    redshift=z,dl=dl)
            return self.zRealiseObject
            #self.ziter=self.redshifts[0] # Need to select z somehow!!
            #self.dataRealisation=lumfuncUtils.calculateH3(cube,self.parameters,\
            #        bins=self.bins,area=self.survey.SURVEY_AREA,family=self.kind,\
            #        redshift=self.ziter)
        else:
            self.dataRealisation=countUtils.calculateI(cube,self.parameters,\
                family=self.kind,bins=self.bins,area=self.survey.SURVEY_AREA,\
                model=self.model)
        return self.dataRealisation

    def transform(self,cube,ndim,nparams):
        #self.params.currentUnit=cube
        #print [c for p,c in enumerate(cube)]
        #for p in self.parameters: self.priorsDict[p]=['U',0.0,1.0]
        #print [cube[i] for i in range(ndim)]
        cube=[self.priors.GeneralPrior(cube[i],self.priorsDict[self.parameters[i]][0],self.priorsDict[self.parameters[i]][1],self.priorsDict[self.parameters[i]][2]) for i in range(ndim)]
        #self.params.currentPhys=cube
        self.setParams(cube)
        return #self.params.currentPhys

    def logprior(self,cube,ndim,nparams):
        #print 'c0',[cube[i] for i in range(ndim)]
        self.currentUnitParams=[cube[i] for i in range(ndim)]
        # NB List comprehension will not modify cube in place
        for i in range(ndim):
            cube[i]=self.priors.GeneralPrior(cube[i],self.priorsDict[self.parameters[i]][0],self.priorsDict[self.parameters[i]][1],self.priorsDict[self.parameters[i]][2])
            #print i,cube[i]
            #if i != 0: cube[i]=1.0
            #print 'c1',[cube[i] for i in range(ndim)]
        self.currentPhysParams=[cube[i] for i in range(ndim)]
        #self.setParams(cube)
        #self.transform(cube,ndim,nparams)
        #print self.currentPhysParams
        return

    def loglike(self,cube,ndim,nparams):
        # Test the break positions if necessary (Si present in params list)
        if not strictly_increasing([cube[i] for i in range(ndim) if \
            (self.parameters[i].startswith('S') and not self.parameters[i].startswith('SL'))]):
            #print [cube[i] for i in range(ndim)]
            print '+',
            return -1.0e99
        else:
            #self.transform(cube,ndim,nparams)
            #print self.currentPhysParams
            if self.doRedshiftSlices:
                if not (cube[self.parameters.index('LMIN')]<cube[self.parameters.index('LSTAR')]<cube[self.parameters.index('LMAX')]):
                    print '+',
                    return -1.0e99
                else:
                    return poissonLhoodMulti2(self.zDataObject,self.realise(cube),\
                                         silent=True)
            elif self.survey.multi:
                return poissonLhoodMulti(self.data,self.realise(cube),\
                                         silent=True,fractions=self.fractions,\
                                         redshifts=self.redshifts)
            else:
               return poissonLhood(self.data,self.realise(cube),silent=True)

    def __str__(self):
        return self.name


#-------------------------------------------------------------------------------
        
    
