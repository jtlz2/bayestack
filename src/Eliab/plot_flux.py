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
rc('legend', fontsize=28)

fig = figure()
offset=[0.6, 0.6, 0.9,  1.,1.25,1.5  ,1.45,2.1,1., 1.05]
for num in range(1,11):
    cat=numpy.loadtxt('cos_data/cos_s%s.txt'%(num))
    cat2=numpy.loadtxt('cos_data/cos_s%s_noise.txt'%(num))
    z_min = zbins[int(num) -1]
    z_max = zbins[int(num)]
    
    fig.add_subplot(4,3,num)
    print 'Flux range/uJy = %f -> %f' % (cat[:,BIN_COL].min(),cat[:,BIN_COL].max())
    #fig = figure()
    binwidth=0.4*SURVEY_NOISE
    #print binwidth
    #hist(cat2[:,BIN_COL], bins=numpy.arange(min(cat[:,BIN_COL]),(20.0*SURVEY_NOISE)+binwidth,binwidth),histtype='step',color='blue')
    n,b,p=hist(cat[:,BIN_COL], bins=numpy.arange(min(cat[:,BIN_COL]),(20.0*SURVEY_NOISE)+binwidth,binwidth),histtype='step',color='blue')
    #yscale('log')
    #xlim(bins[0],20.0*SURVEY_NOISE)
    xlim(-13,21)
    
    #xlabel('S/$\mu$Jy',fontsize=18)
    #ylabel('Number of objects',fontsize=18)
    #coef, unc = curve_fit(myg,b[1:],n,p0=[0.1,1,500])
    #y3 =coef[2]* gaussian(b,coef[0]-5.29755559e-01,coef[1],norm=False)
    y2 = numpy.max(n)*gaussian(b,offset[num-1],SURVEY_NOISE*1.45,norm=False)
    y = numpy.max(n)*gaussian(b,0,SURVEY_NOISE*1.45,norm=False)
    plot(b,y,'r--',linewidth=2,label=r'$\rm{Guassian}$ $\rm{fit}$ $\rm{to}$ $\rm{noise}$')
    plot(b,y2,'b-',linewidth=2,label=r'$\rm{Guassian}$ $\rm{noise}$ $\rm{shifted}$ $\rm{to}$ $\rm{data}$')
    #axvline(1.0*SURVEY_NOISE,color='g',linewidth=2)
    axvline(5.0*3.78,color='g',linewidth=2,linestyle='--',label=r'$\rm{5\sigma}$')
    #text(SURVEY_NOISE,340,'$ \sigma$',rotation=90,color='b',alpha=0.5,fontsize=18)
    #text(5.0*SURVEY_NOISE,340,'$5 \sigma$',rotation=90,color='b',alpha=0.5,fontsize=18)
    #text(-5,20,r'$\rm{%.2f < z < %.2f}$\n $\rm{N = %d}$'%(z_min,z_max,'\n',len(cat)),alpha=0.5,fontsize=18 )
    #tick_params(axis='xaxis',which = 'major', labelsize=25,width =2)
    #tick_params(axis='xaxis',which = 'minor', labelsize=25, width=1)
    #xaxis.set_minor_locator(AutoMinorLocator())
    #if num ==7:
    #    ylabel('Number of objects',fontsize=22)
        
    
    if num==8 or num==9:
    	xticks([-10,-5, 0, 5,10, 15, 20],fontsize=18)
    if num == 10:
    	xticks([-10,-5, 0, 5, 10, 15,20],fontsize=18)
    if num<4:
        xticks( visible=False)
        ylim(0,2550)
        text(-12,2000,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N = %d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
        if num==1:
            yticks([0,600,1200,1800,2400], fontsize=20)
            continue
    elif num<7:
        xticks( visible=False)
        ylim(0,2800)
        text(-12,2100,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N = %d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
        if num==1 or num ==4:
            yticks([0,600,1200,1800,2400], fontsize=20)
            continue
    elif num <10:
        ylim(0,1500)
        text(-12,1200,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N = %d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
        if num==7:
            yticks([0,400,800,1200], fontsize=20)
            continue
    if num==10:
       #if num==1:
       ylim(0,440)
       yticks([0,100,200,300], fontsize=20)
       text(-12,350,r'$\rm{%.2f < z < %.2f}$ %s $\rm{N = %d}$'%(z_min,z_max,'\n',len(cat)),alpha=1,fontsize=22 )
       legend(frameon=False).draggable()
       #else:
       # yticks([1,1e1,1e2])
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

fig.text(0.06,0.5,r'$\rm{Number}$ $\rm{of}$ $\rm{objects}$',fontsize=28, ha='center',va='center', rotation='vertical')   
fig.text(0.5,0.06,r'$\rm{Flux-density (\mu Jy)}$',fontsize=28, ha='center',va='center')
#fig.text(0.5,0.06,r'$\rm{S/(\mu Jy/beam)}$',fontsize=24, ha='center',va='center')
#for a
subplots_adjust(hspace=0,wspace=0)
show()#savefig(BOUT_HISTO)

