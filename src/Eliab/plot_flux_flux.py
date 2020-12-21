#!/usr/bin/env python
from pylab import*
from matplotlib.ticker import AutoMinorLocator
import numpy
from utils import *

bins = [-1056., -750., -500., -400., -350., -300., -250., -200.,\
			-170., -150., -130., -115., -100., -85., -70., -60.,\
			-50., -40., -30., -10., -0.5, 0.5, 10. ,20., 30., 40.,\
			 50., 60., 70, 80., 90., 100., 110., 120., 135., 150.,\
			 165., 180., 195., 210., 225., 240., 260., 280., 300.,\
			 330., 360., 400., 450., 500., 600., 650., 750.]
bins=numpy.array(bins)
zbins =[0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
BIN_COL = -1
SURVEY_NOISE=2.5 # uJy

fig = figure()
offset=[0.6, 0.6, 0.9,  1.,1.25,1.5  ,1.45,2.1,1., 1.05]
for num in range(1,11):
    #cat=numpy.loadtxt('cos_data/cos_flux_s%s.txt'%(num))
    cat=numpy.loadtxt('cos_data/cos_flux_s%s.txt'%(num))
    z_min = zbins[int(num) -1]
    z_max = zbins[int(num)]
    
    fig.add_subplot(4,3,num)
    print 'Flux range/uJy = %f -> %f' % (cat[:,BIN_COL].min(),cat[:,BIN_COL].max())
    #fig = figure()
    plot(cat[:,-1],cat[:,-2],'o')
    plot(cat[:,-1],cat[:,-1],'r')

    xlim(5,2e4)
    ylim(5,2e4)
    xscale('log')
    yscale('log')
    
    
    if num==8 or num==9:
    	xticks([10,100,1000,10000],fontsize=18)
    if num == 10:
    	xticks([10,100,1000,10000],fontsize=18)
    if num<4:
        xticks( visible=False)
        #ylim(0,2550)
        text(11,3000,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N_{VLA} =%d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
        if num==1:
            yticks([10,100,1000,10000], fontsize=20)
            continue
    elif num<7:
        xticks( visible=False)
        #ylim(0,2800)
        text(11,3000,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N_{VLA} =%d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
        if num==1 or num ==4:
            yticks([10,100,1000,10000], fontsize=20)
            continue
    elif num <10:
        #ylim(0,1500)
        text(11,3000,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N_{VLA} =%d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
        if num==7:
            yticks([10,100,1000,10000], fontsize=20)
            xticks( visible=False)
            continue
    if num==10:
       #if num==1:
       #ylim(0,440)
       #yticks([0,10,100,1000,10000], fontsize=20)
       text(11,3000,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N_{VLA} =%d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
       #legend(frameon=False).draggable()
       #else:
       yticks([1e1,1e2,1e3,1e4],fontsize=20)
       continue
    ''' 
      if num==11:
        xlabel(r'$\rm{S/$\mu$Jy',fontsize=22)
        
    if num<10:
        xticks( visible=False)
    if num==10 or num==11:
    	xticks([-1000,0, 1000, 2000],fontsize=16)
    if num == 12:
    	xticks([-1000,0, 1000, 2000, 3000],fontsize=16)
    
    if num==1 or num ==4 or num==7 or num==10:
       #if num==1:
       yticks([1,1e1,1e2], fontsize=16)
       #else:
       # yticks([1,1e1,1e2])
       continue
     '''
       
    yticks(visible=False)  

fig.text(0.06,0.5,r'$\rm{Extracted}$ $\rm{flux-density (\mu Jy)}$',fontsize=28, ha='center',va='center', rotation='vertical')   
fig.text(0.5,0.06,r'$\rm{Catalogue}$ $\rm{flux-density (\mu Jy)}$',fontsize=28, ha='center',va='center')
#for a
subplots_adjust(hspace=0,wspace=0)
show()#savefig(BOUT_HISTO)

