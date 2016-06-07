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
import line_profiler
from utils import *
import matplotlib
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

    if 'labels' in kwargs:
        labels=kwargs['labels']
        parameters=labels # How did this ever work without??
    else:
        labels = ['x', 'y']

    # Furniture fiddle factors
    if 'furniture' in kwargs:
        furniture=kwargs['furniture']
    else:
        furniture={'TRUNCATE_C':False,'TRUNCATE_C_LIMIT':2.0e7,\
               'C_COL':labels.index('C'),'FONTSIZE':4,'ROTATION':60.0,\
               'FIGSIZE':(8.27,11.69), 'DPI':400,\
               'AXIS_LABEL_OFFSET':-0.3,'LOG_BINS':[labels.index('C')]}
    TRUNCATE_C=furniture['TRUNCATE_C']
    TRUNCATE_C_LIMIT=furniture['TRUNCATE_C_LIMIT']
    C_COL=furniture['C_COL']
    FONTSIZE=furniture['FONTSIZE']
    ROTATION=furniture['ROTATION']
    FIGSIZE=furniture['FIGSIZE']
    DPI=furniture['DPI']
    AXIS_LABEL_OFFSET=furniture['AXIS_LABEL_OFFSET']
    log_bins=furniture['LOG_BINS']
    #pylab.gcf().subplots_adjust(left=0.2)

    # !!!! BEWARE THE BINSIZE --- PLOT IS A STRONG FUNCTION OF THIS
    if 'binsize' in kwargs:
        binsize=kwargs['binsize']
    else:
        binsize=50
    print 'Using binsize = %i' % binsize

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

    # Set up the bins
    bins_arr={}
    #log_bins=[labels.index('C')]
    for a in p:
        bins_arr[a]=numpy.linspace(chain[:,a].min(),chain[:,a].max(),binsize)
        if a in log_bins:
            bins_arr[a]=numpy.logspace(numpy.log10(chain[:,a].min()),\
                                       numpy.log10(chain[:,a].max()),binsize)

    # Set up the 2-D panels
    lims={}
    for panel in pairs:
        ipanel+=1
        #bin_spec=(binsize,binsize)
        bin_spec=(bins_arr[panel[0]],bins_arr[panel[1]])
        H,xedges,yedges=numpy.histogram2d(chain[:,panel[0]],chain[:,panel[1]],bins=bin_spec)

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
        assert([N95,N68,N100]==sorted([N95,N68,N100]))
        #Z=reversed(Z); X=reversed(X); Y=reversed(Y)

        levels=[N68,N95,numpy.inf].sort()
        if 'col' in kwargs:
            col=kwargs['col']
        else:
            #col =('#FFFFFF','#a3c0f6','#0057f6') #A pretty blue (pale then dark)
            col =('#a3c0f6','#0057f6')
        #ccol=[x for y, x in sorted(zip(col, [N95,N68,N100]))]
        # Now construct the subplot
        ax[ipanel]=pylab.subplot2grid((nparams,nparams),panel[::-1]) # Reverse quadrant

        #levels=[N68,N95,N100,numpy.inf].sort()
        # Fix levels https://github.com/dfm/corner.py/pull/73
        levelshiftfix=1.0e-4
        levels=[N68-levelshiftfix,N95,N100+levelshiftfix,numpy.inf]
        #for n,l in enumerate(levels[:-2]):
        #    print n,l,col[n]
        #print levels
        levels.sort()
        #print levels
        #print col
        if 'line' in kwargs and kwargs['line']==True:
            CS=pylab.contour(X,Y,Z,levels=levels,colors=col,linewidth=100)
        else:
            CS=pylab.contourf(X,Y,Z,levels=levels,colors=col)


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
            pylab.plot(truth[labels[panel[0]]],truth[labels[panel[1]]],'g+',\
            markersize=20)

        if 'labelDict' in kwargs and kwargs['labelDict'] is not None:
            labelDict=kwargs['labelDict']
        else:
            labelDict=dict((name,name) for name in parameters)

        # Set the axis labels only for left and bottom:
        #print ax[ipanel].get_xlabel(),ax[ipanel].get_ylabel()
        if panel[1] == (nparams-1):
            ax[ipanel].set_xlabel(labelDict[labels[panel[0]]],fontsize=8)
            ax[ipanel].xaxis.set_label_coords(0.5,AXIS_LABEL_OFFSET) # align axis labels
            x_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            ax[ipanel].xaxis.set_major_formatter(x_formatter)
        else:
            ax[ipanel].set_xlabel('')
            ax[ipanel].get_xaxis().set_ticklabels([])

        if panel[0] == 0:
            ax[ipanel].set_ylabel(labelDict[labels[panel[1]]],fontsize=8)
            ax[ipanel].yaxis.set_label_coords(AXIS_LABEL_OFFSET,0.5) # align axis labels
            y_formatter = matplotlib.ticker.ScalarFormatter(useOffset=False)
            ax[ipanel].yaxis.set_major_formatter(y_formatter)
        else:
            ax[ipanel].set_ylabel('')
            ax[ipanel].get_yaxis().set_ticklabels([])

        # Handle log bins
        if panel[0] in log_bins:
            ax[ipanel].set_xscale('log')
            if panel[1] != 0:
                ax[ipanel].get_xaxis().set_ticklabels([])
        elif panel[1] in log_bins:
            ax[ipanel].set_yscale('log')
            if panel[0] != 0:
                ax[ipanel].get_yaxis().set_ticklabels([])

        # Set plot limits
        if autoscale:
            # HACK FOR C ONLY:
            if TRUNCATE_C and panel[0]==C_COL:
                xxlo,xxhi=ax[ipanel].xaxis.get_data_interval()
                if xxhi>TRUNCATE_C_LIMIT:
                    pylab.xlim(xxlo,TRUNCATE_C_LIMIT)
                ax[ipanel].set_xscale('log')
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

        # Save the axis limits (tetris-style)
        lims[panel[0]]=ax[ipanel].xaxis.get_data_interval()
        if panel[1]==nparams-1:
            lims[panel[1]]=ax[ipanel].yaxis.get_data_interval()

        # Some housekeeping
        pylab.xticks(fontsize=FONTSIZE,rotation=ROTATION)
        pylab.yticks(fontsize=FONTSIZE,rotation=0)

    # Set up the 1-D plots on the diagonal
    for iparam in range(nparams):
        bin_spec=bins_arr[iparam]
        J,edges=numpy.histogram(chain[:,iparam],density=True,bins=bin_spec)
        ax1d=pylab.subplot2grid((nparams,nparams),(iparam,iparam))
        pylab.plot(edges[:-1],J,color='k')
        #print iparam,nparams,labels[iparam]

        if 'truth' in kwargs and kwargs['truth'] is not None:
            truth=kwargs['truth']
            pylab.axvline(truth[parameters[iparam]],color='g')

        if iparam == 0:
            ax1d.set_ylabel(labelDict[labels[iparam]],fontsize=8)
            ax1d.yaxis.set_label_coords(AXIS_LABEL_OFFSET,0.5) # align axis labels
        if iparam == (nparams-1):
            ax1d.set_xlabel(labelDict[labels[iparam]],fontsize=8)
            ax1d.xaxis.set_label_coords(0.5,AXIS_LABEL_OFFSET) # align axis labels

        # Set plot limits
        #parameters=['x', 'y', 'S', 'sig', 'Q', 'el', 'em', 'R']
        if autoscale:
            # HACK FOR C ONLY:
            if TRUNCATE_C and iparam==C_COL:
                xxlo,xxhi=ax1d.xaxis.get_data_interval()
                if xxhi>TRUNCATE_C_LIMIT:
                    pylab.xlim(xxlo,TRUNCATE_C_LIMIT)
                #ax1d.set_xscale('log')
                autoscale=True
            xlo,xhi=lims[iparam]
            pylab.xlim(xlo,xhi)
        if not autoscale:
            xlo,xhi=ranges[parameters[iparam]]
            pylab.xlim(xlo,xhi)

        # Handle log bins
        if iparam in log_bins:
            ax1d.set_xscale('log')
        if iparam < (nparams-1):
            ax1d.get_xaxis().set_ticklabels([])
        ax1d.get_yaxis().set_ticklabels([])
        pylab.xticks(fontsize=FONTSIZE,rotation=ROTATION)
        pylab.yticks(fontsize=FONTSIZE)

    #ax1d.set_xscale('linear')

    #axinfo=pylab.subplot2grid((nparams,nparams),(0,nparams-3))
    axinfo=pylab.subplot2grid((nparams,nparams),(0,nparams-nparams%2-1))
    axinfo.get_xaxis().set_visible(False)
    axinfo.get_yaxis().set_visible(False)
    pylab.axis('off')
    pylab.title(title)

    # Plot the truth - this needs to be generalized
    if 'truth' in kwargs and kwargs['truth'] is not None:
        truth=kwargs['truth']
        note=['nparams %i\n truth:' % nparams]
        for k,v in truth.items():
            notelet='%s = %4.2f' % (k,v)
            note.append(notelet)
        pylab.text(0,0,'\n'.join(note))

    if 'outfile' in kwargs:
        outfile=kwargs['outfile']
        pylab.savefig(outfile,figsize=FIGSIZE,dpi=DPI)
        print 'Run: open %s' % outfile
        #pylab.close()
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
