#  Code to plot a contour from an MCMC chain
#  Author: Michelle Knights (2013)

import numpy
import pylab
from scipy import interpolate

#Finds the 95% and 68% confidence intervals, given a 2d histogram of the likelihood
def findconfidence(H):
    H2 = H.ravel()
    H2 = numpy.sort(H2)
    
    #Cut out the very low end
    #H2 = H2[H2>100]

    #Loop through this flattened array until we find the value in the bin which contains 95% of the points
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

#Given a chain, labels and a list of which parameters to plot, plots the contours
# Arguments:
# chain=an array of the chain (not using weights, i.e. each row counts only once)
# p= a list of integers: the two parameters you want to plot (refers to two columns in the chain)
#kwargs: labels= the labels of the parameters (list of strings)
#               col=a tuple of the two colours for the contour plot
#               line=boolean whether or not to just do a line contour plot
def contour(chain,p,**kwargs):
    binsize=50
    H, xedges, yedges = numpy.histogram2d(chain[:,p[0]],chain[:,p[1]], bins=(binsize,binsize))
    
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

    if 'line' in kwargs and kwargs['line']==True:
        pylab.contour(X, Y,Z,levels=[N95,N68,N100],colors=col, linewidth=100)
    else:
        pylab.contourf(X, Y,Z,levels=[N95,N68,N100],colors=col)
    pylab.xlabel(labels[p[0]],fontsize=22)
    pylab.ylabel(labels[p[1]],fontsize=22)
    pylab.show()

#Testing all functionality
#c=pylab.loadtxt('chain_2d_banana.txt')
#contour(c,[0,1], labels=['1', '2'], col=('#3bf940','#059a09'),line=True)
