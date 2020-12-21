#!/usr/bin/env python

"""
This is reconstruct_lf.py
June 20 
Eliab based on Jon's reconstruct.py

Usage:
Works for LFs (Sch,DPL,DDPL)
./reconstruct.py CHAINS_DIR

"""

import os,sys
import importlib
import numpy
from scipy import stats
import pylab as plt
import pymultinest
from bayestackClasses import countModel
from utils import sqDeg2sr,calculate_confidence,peak_confidence,fetchStats,calculate_confidence2
from lumfuncUtils import get_sbins,get_z,get_dl,get_sfr_q, get_sfrd_z, get_q,sfrd_Behroozi, sfrd_Madau, LF, get_Lbins, doublepowerlaw

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

    # Set up the experiment
    #bins=[12.679592, 31.84969509, 80.00281696, 200.9579904, 504.78364938, 1010, 3267.95919972]# 0.4 bins
    
    bins=numpy.logspace(0,3.01,10)
    bins = list(bins)
    bins.append(2e3)
    bins2 = numpy.arange(18.,27.,0.2)
    bins1 = numpy.arange(18.2,27.2,0.2)
    print redshifts
    nbins = len(bins)
    z_m = redshifts[0]
    dl = get_dl(z_m) 
    sbin1 = get_sbins(10**bins1,z_m,dl)*1e6
    sbin2 = get_sbins(10**bins2,z_m,dl)*1e6
    #L = 10**bins1*(1.4/3.)**(.7)
    #L =get_Lbins(sbin1,z_m,dl)*(1.4/3.)**(-.7)
    expt=countModel(modelFamily,nlaws,settingsf,[dataset],floatNoise,doRedshiftSlices=True,mybins=sbin2)
    f='%spost_equal_weights.dat' % outstem
    f=os.path.join(outdir,f)    
    print 'os.pathf',f
    chain=['a','b','c','d','e','f','g','h','i','j','k']
     # Load equally-weighted posterior samples
    x=numpy.genfromtxt(f)
    nsamp=x.shape[0]
    ncols=x.shape[1] # The fifth [seventh] column is the posterior value
     # There must be a better way, but:
     #ncols = 14
    z=numpy.zeros((nsamp,ncols-1+expt.nbins))
    z[:,:-(expt.nbins-1)]=x
     # Shift posterior values to end
    z[:,-1]=z[:,ncols-1] # Copy...
    z[:,ncols-1]=0.0     # ...and blank
     
    print 'outstem',outstem

     # Fetch best-fit parameters and calculate best-fit line
    ana=pymultinest.analyse.Analyzer(ncols-1,\
        outputfiles_basename=os.path.join(outdir,outstem))
    drawmap=ana.get_best_fit()['parameters']
        
    if True:
        print '--> Calculating *ML* reconstruction'
        drawmap=drawmap

    ymap=expt.evaluate(drawmap)#,expt.binsMedian)
    #print ymap

    evaluations=[]   
    evaluate=True 
    for num in range(1,12):
     print num
     
     z_min,z_max, z_m = get_z(num,False)
      
    # for k in range(len(phi_sf)):
    #    print sbin1[k],bins1[k],numpy.power(10,bins1[k]),numpy.log10(L_1[k]),numpy.log10(L_2[k]), numpy.log10(phi_sf[k]), numpy.log10(phi_agn[k]),numpy.log10(ymap[k]),numpy.log10(ymap2[k]),numpy.log10(ymap_ev[num-1][k])
     #sys.ek    
     #print ymap
                     
     for isamp in xrange(nsamp):
        if evaluate:
            evaluations.append(expt.evaluate(z[isamp,:]))
            
        z[isamp,ncols-1:]=evaluations[isamp][num-1]
                          #[LF(S,z_m,0,0,dl,z[isamp,:],\
                          #expt.parameters,area=expt.survey.SURVEY_AREA,\
                          #family=modelFamily) for S in expt.binsMedian/1.0e6]
        #print             [LF(S,z_m,0,0,dl,z[isamp,:],\
        #                  expt.parameters,area=expt.survey.SURVEY_AREA,\
        #                  family=modelFamily) for S in expt.binsMedian/1.0e6]
        #sys.exit()
        
     #sys.exit()
     evaluate=False
     z[numpy.where(z==0.0)]='NaN'
    
     s=numpy.zeros((len(expt.binsMedian),8))
     s[:,0]=expt.binsMedian
     recons=z[:,ncols-1:]

    #print '# ibin flux fit low high dlower dupper skew kurtosis'
     print '%6s %12s %12s %12s %12s %12s %12s %12s %12s'%('bin', 'lower phi_rms','lower phi', 'l_err','ymap','u_err', 'upper phi','u_phi_rms', '2*rms')
     for ibin in xrange(len(expt.binsMedian)):
        x = recons[:,ibin]
        # Remove NaNs from stats vectors
        # http://stackoverflow.com/questions/11620914/removing-nan-values-from-an-array
        x = x[~numpy.isnan(x)]
        #ss=stats.bayes_mvs(x,alpha=0.68)[0]
        #x*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print bins1[ibin]
        #plt.hist(numpy.log10(x),bins=15)
        #plt.xlabel(r'$\rm{log_{10}[\rho_m(Mpc^{-3} mag^{-1})]}$',fontsize = 20)
        #plt.show()
        tt=peak_confidence(x,bins=10)

        try:
            ss=numpy.zeros(3)
            ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence(x,alpha=0.95,ret_all=True)
            #ss[0],dlow,dhigh,ss[1],ss[2]=calculate_confidence2(x,alpha=0.68,ret_all=True,\
            #                                                   value_central=ymap[ibin],\
            #                                                   truncate_edges=True)
        except:
            ss=numpy.nan*numpy.ones(3)
            print "didn't work for ", x
            continue
        #sys.exit()
        #tt=peak_confidence(x,bins=10)
        #ss*=numpy.power(s[ibin,0]/1.0e6,2.5)
        #print ss[0],tt
#        s[ibin,1]=ss[0]  # median
#        s[ibin,1]=tt     # peak
        s[ibin,1]=ymap[num-1][ibin] # MAP
        #print ymap
        s[ibin,2]=ss[1]  # lower
        s[ibin,3]=ss[2]  # upper
        s[ibin,4]= 2*numpy.std(x)# (s[ibin,3]- s[ibin,2])/(2*numpy.std(x))# skewness
        s[ibin,5]= ss[0]
        
        #print '%6.1f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f %12.3f'%(bins1[ibin], numpy.log10((ymap[ibin] - s[ibin,4])), numpy.log10(ss[1]),numpy.log10(s[ibin,2]),numpy.log10(ymap[ibin]),numpy.log10(s[ibin,3]) ,numpy.log10(ss[2]), numpy.log10(s[ibin,4]+ ymap[ibin]), numpy.log10(s[ibin,4]) )
        #print ibin,s[ibin,0],s[ibin,1],dlow,dhigh,ss[1],ss[2],s[ibin,4],s[ibin,5]#,stats.skewtest(x)
        #sys.exit()

    # ...and output to file
    #sys.exit()
     rstatsf='recon_stats_%s.txt'%(chain[num-1])
     rstatsf=os.path.join(outdir,rstatsf)
     hdr='# median_flux_uJy dnds_2p5_Jy1p5srm1 delta_dnds_2p5_lower_Jy1p5srm1 delta_dnds_2p5_upper_Jy1p5srm1 skewness kurtosis'
     fid = open(rstatsf,'w')
     print hdr
     print s

     print '-> Writing stats (i.e. dnds2p5+/- v. S/uJy) to %s' % rstatsf
     fid.write('%s\n'%hdr)
     numpy.savetxt(fid,s)
     fid.close()
     
     print 'Finished.'
    


    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)
