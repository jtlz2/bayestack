import os
import importlib
import numpy
import itertools
from scipy import integrate
from scipy.special import erf
from priors import Priors
import countUtils
from utils import sqDeg2sr,beamFac,sqrtTwo,strictly_increasing,poissonLhood
import cosmolopy

#-------------------------------------------------------------------------------

class surveySetup(object):

    """
    Set up the survey parameters (assumed fixed)
    survey=surveySetup(whichSurvey)
    """

    def __init__(self,whichSurvey):
        self.whichSurvey=whichSurvey
        if whichSurvey == 'video':
            self.datafile=os.path.join(whichSurvey,'all_test_41_150120a.txt')
            self.HALO_MASK=11436315.0/(19354.0*19354.0)
            self.SURVEY_AREA=1.0 *(1.0-self.HALO_MASK)# sq.deg. [Boris -> 0.97 sq. deg.]
            self.SURVEY_NOISE=16.2 # uJy [median; mean=16.3]
            self.radioSynthBeamFWHM=4.0 # pixels/upsamplingFactor
            self.radioSynthOmegaSr=sqDeg2sr*beamFac*(self.radioSynthBeamFWHM/3600.0)**2

#-------------------------------------------------------------------------------

class binSetup(object):
    """
    binScheme=binSetup(whichBins)
    """
    
    def __init__(self,whichBins):
        self.whichBins=whichBins
        bins=numpy.array([-108.0,-80.0,-50.0,-30.0,-20.0,-10.0,-5.0,-2.5,-1.0,-0.5,0.0,0.5,1.0,2.5,5.0,7.5,10.0,15.0,20.0,25.0,30.0,40.0,50.0,65.0,85.0])

        self.bins=bins
        nbins=len(bins)-1
        self.nbins=nbins

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
            self.func=cosmolopy.utils.PiecewisePowerlaw(params['breaks'],\
                        params['slopes'],externalval=0.0,coefficients=None)
        elif self.family=='poly':
            self.func=numpy.poly1d(params['coeffs'])

        return params['amp']*self.func(vals)

#-------------------------------------------------------------------------------

class countModel(object):

    """
    usage:
    counts=countModel(kind,order,settingsf,whichSurvey,whichBins)
    
    """

    def __init__(self,kind,order,settingsf,whichSurvey,whichBins,floatNoise):
        # Import settings
        set_module=importlib.import_module(settingsf)
        globals().update(set_module.__dict__)

        # Set up model
        self.kind=kind
        self.name=kind
        self.order=order
        self.nlaws=order
        self.model=model(self.kind)
        self.floatNoise=floatNoise

        # Set up parameters for this model
        self.paramsAvail=\
                    {'breaks':['S%i'%ic for ic in xrange(self.nlaws+1)],\
                     'slopes':['a%i'%ic for ic in xrange(self.nlaws)],\
                     'coeffs':['p%i'%ic for ic in xrange(self.nlaws)],\
                     'amp':['C'],'extra':['noise']}
        familyMap={'ppl':['breaks','slopes','amp','extra'],'poly':['coeffs','extra']}
        self.paramsStruct=\
          [self.paramsAvail[p] for p in self.paramsAvail if p in familyMap[kind]]
        # --> This defines the order of the parameters:
        self.parameters=list(itertools.chain(*self.paramsStruct))
        self.nparams=len(self.parameters)
        print self.nparams,self.parameters
        self.currentPhysParams=-99.0*numpy.ones(self.nparams)

        # Set up priors
        self.priors=Priors()
        self.priorsDict=self.parsePriors(self.parameters,self.floatNoise)

        # Set up data and bins
        self.survey=surveySetup(whichSurvey)
        self.binScheme=binSetup(whichBins)
        self.bins=self.binScheme.bins
        self.nbins=self.binScheme.nbins
        self.data,self.bins=self.loadData(self.survey.datafile)

        return

    def parsePriors(self,parameters,floatNoise):
        priorsDict={}
        for p in parameters:
            if p[0]=='C': priorsDict[p]=['LOG',C_MIN,C_MAX]
            elif p[0]=='S': priorsDict[p]=['U',SMIN_MIN,SMAX_MAX]
            elif p[0]=='a': priorsDict[p]=['U',ALPHA_MIN,ALPHA_MAX]
            elif p[0]=='n':
                if floatNoise:
                    priorsDict[p]=['U',0.5*SURVEY_NOISE,2.0*SURVEY_NOISE]
                else:
                    priorsDict[p]=['DELTA',SURVEY_NOISE,SURVEY_NOISE]
            if p=='S0': priorsDict[p]=['U',SMIN_MIN,SMIN_MAX]
            if p=='S1': priorsDict[p]=['U',SMAX_MIN,SMAX_MAX] # --> generalize
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

    def evaluate(self,Ss,params):
        return self.model.eval(Ss,self.currentParams)

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

    def dn_by_ds(self,return_all):
        return countUtils.calculateDnByDs(self.bins,self.dataRealisation,\
                                return_all=return_all)

    def erfs(self,S,Sbin,sigma):
        """vectorized"""
        Sbinlow,Sbinhigh=Sbin
        e=0.5*(erf((S-Sbinlow)/(sqrtTwo*sigma)) - erf((S-Sbinhigh)/(sqrtTwo*sigma)))
        return e

    def realise(self,cube):
        self.dataRealisation=countUtils.calculateI3(cube,self.parameters,\
                family=self.kind,bins=self.bins,area=self.survey.SURVEY_AREA)
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
        #print 'c1',[cube[i] for i in range(ndim)]
        self.currentPhysParams=[cube[i] for i in range(ndim)]
        #self.setParams(cube)
        #self.transform(cube,ndim,nparams)
        #print self.currentPhysParams
        return

    def loglike(self,cube,ndim,nparams):
        # Test the break positions if necessary
        if self.kind=='ppl' and not strictly_increasing([cube[i] for i in range(ndim) if self.parameters[i][0]=='S']):
            print '+',
            return -1.0e99
        else:
            #self.transform(cube,ndim,nparams)
            #print self.currentPhysParams
            return poissonLhood(self.data,self.realise(cube))

    def __str__(self):
        return self.name


#-------------------------------------------------------------------------------
        
    
