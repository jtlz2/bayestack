#!/usr/bin/env python
 
"""
Code to plot a contour from an MCMC chain
Author: Michelle Knights (2013)
Modified: Jonathan Zwart (12 August 2013)
"""
 
import sys,os
import numpy
import pylab
from scipy import interpolate
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#from lumfunc import *
import line_profiler
from utils import *
#from settings import *
 
from matplotlib.path import Path
 
 
#-------------------------------------------------------------------------------
 
 
def findconfidence(H):
    """
    Finds the 95% and 68% confidence intervals, given a 2d histogram
    of the likelihood
    """
 
 
    H2 = H.ravel()
    H2 = numpy.sort(H2)
     
    #Cut out the very low end
    #H2 = H2[H2>100]
 
    #Loop through this flattened array until we find the value in the
    #bin which contains 95% of the points
    tot = sum(H2)
    tot95=0
    tot68=0
 
    #Changed this to 68% and 30% C.I
    for i in range(len(H2)):
        tot95 += H2[i]
        if tot95 >= 0.05*tot:
            N95 = H2[i]
            #print i
            break
 
    for i in range(len(H2)):
        tot68 += H2[i]
        if tot68>=0.32*tot:
            N68 = H2[i]
            break  
    return max(H2),N95,N68
 
#-------------------------------------------------------------------------------
 
def contour(chain,p,**kwargs):
    """
    Original alias for contourSingle
    """
    return contourSingle(chain,p,**kwargs)
 
#-------------------------------------------------------------------------------
 
def contourSingle(chain,p,**kwargs):
    """
    #Given a chain, labels and a list of which parameters to plot, plots the contours
    # Arguments:
    # chain=an array of the chain (not using weights, i.e. each row counts only once)
    # p= a list of integers: the two parameters you want to plot (refers to two columns in the chain)
    #kwargs:        labels= the labels of the parameters (list of strings)
    #               col=a tuple of the two colours for the contour plot
    #               line=boolean whether or not to just do a line contour plot
    #               outfile='outf.png'
    """
 
    # !!!! BEWARE THE BINSIZE --- PLOT IS A STRONG FUNCTION OF THIS
    binsize=50
    H,xedges,yedges=numpy.histogram2d(chain[:,p[0]],chain[:,p[1]],bins=(binsize,binsize))
     
    x=[]
    y=[]
    z=[]
    for i in range(len(xedges[:-1])):
        for j in range(len(yedges[:-1])):
            x.append(xedges[:-1][i])
            y.append(yedges[:-1][j])
            z.append(H[i, j])
 
    SMOOTH=False
    if SMOOTH:
        sz=50
        smth=80e6
        spl = interpolate.bisplrep(x, y, z,  s=smth)
        X = numpy.linspace(min(xedges[:-1]), max(xedges[:-1]), sz)
        Y = numpy.linspace(min(yedges[:-1]), max(yedges[:-1]), sz)
        Z = interpolate.bisplev(X, Y, spl)
    else:
        X=xedges[:-1]
        Y=yedges[:-1]
        Z=H
 
    #I think this is the weird thing I have to do to make the contours work properly
    X1=numpy.zeros([len(X), len(X)])
    Y1=numpy.zeros([len(X), len(X)])
    for i in range(len(X)):
        X1[ :, i]=X
        Y1[i, :]=Y
    X=X1
    Y=Y1
     
    N100,N95,N68 = findconfidence(Z)
 
    if 'col' in kwargs:
        col=kwargs['col']
    else:
        col =('#a3c0f6','#0057f6') #A pretty blue
 
    if 'labels' in kwargs:
        labels=kwargs['labels']
    else:
        labels = ['x', 'y']
 
    pylab.clf()
 
    if 'line' in kwargs and kwargs['line']==True:
        pylab.contour(X, Y,Z,levels=[N95,N68,N100],colors=col, linewidth=100)
    else:
        pylab.contourf(X, Y,Z,levels=[N95,N68,N100],colors=col)
    pylab.xlabel(labels[p[0]],fontsize=22)
    pylab.ylabel(labels[p[1]],fontsize=22)
 
     
    if 'outfile' in kwargs:
        outfile=kwargs['outfile']
        pylab.savefig(outfile)
        #pylab.close()
    else:
        pylab.show()
 
    return
 
#-------------------------------------------------------------------------------
 
 
def contourTri(chain,**kwargs):
    """
    #Given a chain, labels and a list of which parameters to plot, plots the contours
    # Arguments:
    # chain=an array of the chain (not using weights, i.e. each row counts only once)
    # p= a list of integers: the two parameters you want to plot (refers to two columns in the chain)
    #kwargs:        labels= the labels of the parameters (list of strings)
    #               col=a tuple of the two colours for the contour plot
    #               line=boolean whether or not to just do a line contour plot
    #               outfile='triangle.png'
    #               binsize=50
    #               reconstruct=boolean whether or not to plot reconstruction
    #               autoscale=boolean whether or not to autoscale axes
    #               ranges=dictionary of plot range lists, labelled by
    #                      parameter name, e.g. {'A':[0.0,1.0],etc.}
    #               title=outdir
    p is now ignored
    """
 
    # Collate the contour-region info
    bundle=chain
 
 
    TRUNCATE_C=False
    FONTSIZE=15; ROTATION=60.0
    FIGSIZE=(8.27,11.69); DPI=600
    #FIGSIZE=(8.27,11.69); DPI=400
    AXIS_LABEL_OFFSET=-0.5    #-0.7
    
    X_LABEL_OFFSET=-0.35#-0.35    #-0.8
    Y_LABEL_OFFSET=-0.4#-0.4    #-0.6
 
    # !!!! BEWARE THE BINSIZE --- PLOT IS A STRONG FUNCTION OF THIS
    if 'binsize' in kwargs:
        binsize=kwargs['binsize']
    else:
        binsize=50
    print 'Using binsize = %i' % binsize
 
    if 'labels' in kwargs:
        labels=kwargs['labels']
        parameters=labels # How did this ever work without??
    else:
        labels = ['x', 'y']
 
    if 'ranges' in kwargs:
        ranges=kwargs['ranges']
    else:
        ranges=None
 
    if 'title' in kwargs:
        title=kwargs['title']
    else:
        title=''
 
    if 'autoscale' in kwargs:
        autoscale=kwargs['autoscale']
    else:
        autoscale=True
         
    p = range(len(labels))
    pairs = trianglePairs(p)
    nparams = len(p)
 
    # Start setting up the plot
    ipanel=0; ax={}
    pylab.clf()
    for panel in pairs:
        ipanel+=1       
        H,xedges,yedges=numpy.histogram2d(chain[:,panel[0]],chain[:,panel[1]],bins=(binsize,binsize))
 
        x=[]
        y=[]
        z=[]
        for i in range(len(xedges[:-1])):
            for j in range(len(yedges[:-1])):
                x.append(xedges[:-1][i])
                y.append(yedges[:-1][j])
                z.append(H[i, j])
 
        SMOOTH=False
        if SMOOTH:
            sz=50
            smth=80e6
            spl = interpolate.bisplrep(x, y, z,  s=smth)
            X = numpy.linspace(min(xedges[:-1]), max(xedges[:-1]), sz)
            Y = numpy.linspace(min(yedges[:-1]), max(yedges[:-1]), sz)
            Z = interpolate.bisplev(X, Y, spl)
        else:
            X=xedges[:-1]
            Y=yedges[:-1]
            Z=H
     
        #I think this is the weird thing I have to do to make the contours work properly
        X1=numpy.zeros([len(X), len(X)])
        Y1=numpy.zeros([len(X), len(X)])
        for i in range(len(X)):
            X1[ :, i]=X
            Y1[i, :]=Y
        X=X1
        Y=Y1
     
        N100,N95,N68 = findconfidence(Z)
 
        if 'col' in kwargs:
            col=kwargs['col']
        else:
            col =('#a3c0f6','#0057f6') #A pretty blue
 
        # Now construct the subplot
        ax[ipanel]=pylab.subplot2grid((nparams,nparams),panel[::-1]) # Reverse quadrant
 
        if 'line' in kwargs and kwargs['line']==True:
            CS=pylab.contour(X, Y,Z,levels=[N95,N68,N100],colors=col, linewidth=100)
        else:
            CS=pylab.contourf(X, Y,Z,levels=[N95,N68,N100],colors=col)
 
 
        # Identify points lying within 68 percent contour region
        #print ipanel,bundle.shape,chain.shape
        #v=CS.collections[0].get_paths()[0].vertices
        #w=CS.collections[1].get_paths()[0].vertices
        #print v[:,0]
        #print v[:,1]
        #b=bundle[:,[panel[0],panel[1]]]
        #mask=Path(v).contains_points(b)
        #mask2=Path(w).contains_points(b)
        #print panel[0],panel[1],b[:,0].size,b[:,0][mask].size,b[:,0][mask2].size,labels[panel[0]],labels[panel[1]]
 
             
        if 'truth' in kwargs and kwargs['truth'] is not None:
            truth=kwargs['truth']
            #pylab.plot(truth[labels[panel[0]]],truth[labels[panel[1]]],'r+',\
            #markersize=20)
 
        if 'labelDict' in kwargs and kwargs['labelDict'] is not None:
            labelDict=kwargs['labelDict']
        else:
            labelDict=dict((name,name) for name in parameters)
 
        # Set the axis labels only for left and bottom:
        #print ax[ipanel].get_xlabel(),ax[ipanel].get_ylabel()
        if panel[1] == (nparams-1):
            ax[ipanel].set_xlabel(labelDict[labels[panel[0]]],fontsize=18)
            ax[ipanel].xaxis.set_label_coords(0.5,X_LABEL_OFFSET) # align axis labels
            ax[ipanel].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax[ipanel].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
 
        else:
            ax[ipanel].set_xlabel('')
            ax[ipanel].get_xaxis().set_ticklabels([])
        if panel[0] == 0:
            ax[ipanel].set_ylabel(labelDict[labels[panel[1]]],fontsize=18)
            ax[ipanel].yaxis.set_label_coords(Y_LABEL_OFFSET,0.5) # align axis labels
            #ax[ipanel].xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            ax[ipanel].yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        else:
            ax[ipanel].set_ylabel('')
            ax[ipanel].get_yaxis().set_ticklabels([])
 
        # Set plot limits
        if autoscale:
            # HACK FOR C ONLY:
            if TRUNCATE_C and panel[0]==0:
                xxlo,xxhi=ax[ipanel].xaxis.get_data_interval()
                if xxhi>1.0e4:
                    pylab.xlim(xxlo,1.0e4)
                #ax[ipanel].set_xscale('log')
                autoscale=True
                #locs,labels = plt.xticks()
                #plt.xticks(locs, map(lambda x: "%g" % x, locs*1.0e5))
        else:
            xlo=ranges[labels[panel[0]]][0]
            xhi=ranges[labels[panel[0]]][1]
            ylo=ranges[labels[panel[1]]][0]
            yhi=ranges[labels[panel[1]]][1]
            pylab.xlim(xlo,xhi)
            pylab.ylim(ylo,yhi)
 
        # Some housekeeping
        pylab.xticks(fontsize=FONTSIZE,rotation=ROTATION)
        pylab.yticks(fontsize=FONTSIZE,rotation=0)
        pylab.locator_params(nbins=4)
 
    # Set up the 1-D plots on the diagonal
    for iparam in range(nparams):
        #        b=numpy.histogram(R,bins=bins)
        J,edges=numpy.histogram(chain[:,iparam],density=True,bins=binsize)
        ax1d=pylab.subplot2grid((nparams,nparams),(iparam,iparam))
        pylab.plot(edges[:-1],J,color='k')
        #print iparam,nparams,labels[iparam]
 
        if 'truth' in kwargs and kwargs['truth'] is not None:
            truth=kwargs['truth']
            #pylab.axvline(truth[parameters[iparam]],color='g')
 
        if iparam == 0:
            ax1d.set_ylabel(labelDict[labels[iparam]],fontsize=18)
            ax1d.yaxis.set_label_coords(Y_LABEL_OFFSET,0.5) # align axis labels
        if iparam == (nparams-1):
            ax1d.set_xlabel(labelDict[labels[iparam]],fontsize=18)
            ax1d.xaxis.set_label_coords(0.5,X_LABEL_OFFSET) # align axis labels
 
        # Set plot limits
        #parameters=['x', 'y', 'S', 'sig', 'Q', 'el', 'em', 'R']
        if autoscale:
            # HACK FOR C ONLY:
            if TRUNCATE_C and iparam==0:
                xxlo,xxhi=ax1d.xaxis.get_data_interval()
                if xxhi>1.0e4:
                    pylab.xlim(xxlo,1.0e4)
                #ax1d.set_xscale('log')
                autoscale=True
        if not autoscale:
            xlo,xhi=ranges[parameters[iparam]]
            pylab.xlim(xlo,xhi)
            if TRUNCATE_C and iparam==0:
                xxlo,xxhi=ax1d.xaxis.get_data_interval()
                if xxhi>1.0e4:
                    pylab.xlim(xxlo,1.0e4)
 
        if iparam < (nparams-1):
            ax1d.get_xaxis().set_ticklabels([])
        ax1d.get_yaxis().set_ticklabels([])
        pylab.xticks(fontsize=FONTSIZE,rotation=ROTATION)
        pylab.yticks(fontsize=FONTSIZE)
        pylab.locator_params(nbins=4)
        #if iparam == 0: ax1d.set_xscale('log')
    #ax1d.set_xscale('linear')
 
    #axinfo=pylab.subplot2grid((nparams,nparams),(0,nparams-3))
    axinfo=pylab.subplot2grid((nparams,nparams),(0,nparams-nparams%2-1))
    axinfo.get_xaxis().set_visible(False)
    axinfo.get_yaxis().set_visible(False)
    pylab.axis('off')
    pylab.title(title)
 
    # Plot the truth - this needs to be generalized for non-lumfunc
    if 'truth' in kwargs and kwargs['truth'] is not None:
        truth=kwargs['truth']
        note=['nparams %i\n truth:' % nparams]
        for k,v in truth.items():
            notelet='%s = %4.2e' % (k,v)
            note.append(notelet)
        #pylab.text(-1,-1,'\n'.join(note))
 
        if 'reconstruct' in kwargs:
            reconstruct=kwargs['reconstruct']
            axrecon=pylab.subplot2grid((nparams,nparams),(0,nparams-2),\
                                       rowspan=2,colspan=2)
            axrecon.set_xscale('log')
            axrecon.set_yscale('log')
            pylab.xticks(fontsize=FONTSIZE,rotation=60)
            pylab.yticks(fontsize=FONTSIZE)
            pylab.locator_params(nbins=5)
            median_bins=medianArray(reconstruct[0])
            dnds=calculateDnByDs(median_bins,reconstruct[1])
            dndsN=calculateDnByDs(median_bins,ksNoisy)
            print median_bins
            print dnds
            print 'truth items \n'
            print truth.items()
            recon=numpy.zeros(numpy.shape(median_bins))
            post=numpy.zeros(numpy.shape(median_bins))
            print '# i Smedian ks dnds dndsS2.5 NRecon dndsRecon dndsS2.5Recon log10dnds log10dndsR diffR dndsN'
            if nparams == 4:
                (C,alpha,Smin,Smax)\
                  =(truth['C'],truth['alpha'],truth['Smin'],truth['Smax'])
                area=10.0 # Hack
                # Reconstruct powerLaw points given truth
                for i in range(len(median_bins)):
                    recon[i]=powerLawFuncS(median_bins[i],\
                                                   C,alpha,Smin,Smax,area)
                    post[i]=powerLawFuncS(median_bins[i],\
                                                  9.8,-0.63,0.04,14.1,area)
                #recon *= lumfunc.ksRaw
                #dndsR=lumfunc.calculateDnByDs(median_bins,recon)
                # **** XXXX Where does the 1000 come from!? :(( XXXX
                dndsR=recon*1000.0
                dndsP=post*1000.0
                # cols: i Smedian ks dnds dndsS2.5 NRecon dndsRecon
                # dndsS2.5Recon log10dnds log10dndsR diffR dndsN
                for i in range(len(median_bins)):
                    print '%i %f %i %i %i %i %i %i %f %f %i %i' % (i,median_bins[i],\
                                                  reconstruct[-1][i],dnds[i],\
                      dnds[i]*median_bins[i]**2.5,recon[i],dndsR[i],\
                      dndsR[i]*median_bins[i]**2.5,numpy.log10(dnds[i]),\
                        numpy.log10(dndsR[i]),int(dndsR[i]-dnds[i]),dndsN[i])
 
                      #print recon
            pylab.xlim(reconstruct[0][0],reconstruct[0][-1])
            #pylab.ylim(1.0e2,1.0e8)
 
            #pylab.plot(median_bins,dnds*numpy.power(median_bins,2.5)*lumfunc.sqDeg2sr,'+')
            power=2.5
            pylab.plot(median_bins,dnds*sqDeg2sr*numpy.power(median_bins,power),'+')
            pylab.plot(median_bins,dndsR*sqDeg2sr*numpy.power(median_bins,power),'-')
            pylab.plot(median_bins,dndsN*sqDeg2sr*numpy.power(median_bins,power),'+')
            pylab.plot(median_bins,dndsP*sqDeg2sr*numpy.power(median_bins,power),'-')
            #pylab.plot(dnds,dndsR*numpy.power(median_bins,1.0))
            #b=lumfunc.simtable(lumfunc.bins,a=-1.5,seed=1234,noise=10.0,dump=False)
 
    if 'outfile' in kwargs:
        outfile=kwargs['outfile']
        pylab.subplots_adjust(wspace=0, hspace=0)
        #pylab.show()
        pylab.savefig(outfile,figsize=FIGSIZE,dpi=DPI)
        print 'Run: open %s' % outfile
        pylab.close()
    else:
        pylab.show()
 
    return bundle
 
#-------------------------------------------------------------------------------
 
 
def trianglePairs(inlist):
    """
    """
    pairs=[]
    for i in inlist:
        for j in inlist:
            if j > i:
                pairs.append((i,j))
 
    return pairs
     
#-------------------------------------------------------------------------------
 
if __name__ == '__main___':
    """
    """
 
 
    #parameters=lumfunc.parameters['C','alpha','Smin','Smax']
 
#Testing all functionality
#c=pylab.loadtxt('chain_2d_banana.txt')
#contour(c,[0,1], labels=['1', '2'], col=('#3bf940','#059a09'),line=True)
 
    # Run as e.g.
    contour_plot.contourTri(pylab.loadtxt('chains-4-all-10deg-130812/1-post_equal_weights.dat'),line=True,outfile='chains-4-all-10deg-130812/test.png',col=('red','blue'),labels=lumfunc.parameters,ranges=lumfunc.plotRanges,truth=lumfunc.plotTruth,reconstruct=(lumfunc.medianArray(lumfunc.bins),lumfunc.ksNoisy),autoscale=False,title='title')
 
    #contour_plot.contourTri(pylab.loadtxt('chains-3-fixSmin-10deg-130812/1-post_equal_weights.dat'),line=True,outfile='test.png',col=('red','blue'),labels=lumfunc.parameters)
 
    #import pymultinest
    #a=pymultinest.Analyzer(len(lumfunc.parameters),'chains-4-all-mm-10deg-130815/1-')
    #s=a.get_stats()
    #print s
     
    sys.exit(0)
