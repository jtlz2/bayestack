#!/usr/bin/env python
from pylab import*
from matplotlib.ticker import AutoMinorLocator
import numpy
from lumfuncUtils import get_dl, get_Vmax,get_dsdl,get_z
from utils import *
from cosmocalc import cosmocalc

Ho=70.0
wm=0.3

zbins =[0.1, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.2, 4.0]
t= [cosmocalc(zi,H0=70,WM = 0.3)['zage_Gyr'] for zi in zbins]
BIN_COL = -1
SURVEY_NOISE=2.5 # uJy
rc('legend', fontsize=23)

def sfr(m,z,So=2.8,a1=0.23,a2=0.13,mo=10.8):
    return So - a1*z-numpy.log10(1+(10**(mo-a2*z)/10**(m)))
    
def phi(m,m_star,phi1,phi2,a,b):
    m_m=10**(m-m_star)
    #print m_star, phi1,phi2,a,b
    sk=numpy.log(10)* (phi1*(m_m)**(a+1) + phi2*(m_m)**(1+b))
    mf=numpy.exp(-m_m)*sk
    #print numpy.exp(-m_m)
    return mf

#5 0.28<z<0.36, 7 0.45<z<0.56 9 0.68<z<0.82 10 0.82<z<1 11 1<z<1.2 13 1.45<z<1.75
#14 1.75<z<2.2 15 2.2<z<2.6 16 2.6z3.25 17 3.25z3.75
w_18=[[10.88 , -0.55 , -1.536, -2.555, -3.4  , 8.388, 1.592],\
      [11.065, -0.487, -1.423, -3.632, -3.241, 8.115, 0.904],\
      [10.881, -0.858, -1.838, -2.84 , -3.977, 8.204, 2.495],\
      [11.06 , -0.986, -1.591, -3.321, -3.532, 8.138, 3.799],\
      [10.862, -0.522, -1.488, -3.07 , -3.409, 8.025, 2.564],\
      [10.896, -0.36 , -1.653, -3.295, -3.785, 7.891, 6.656],\
      [11.041, -0.351, -1.589, -3.478, -3.648, 7.98 , 0.342],\
      [11.075, -0.336, -1.58 , -4.039, -3.62 , 7.846, 1.008],\
      [11.471, -0.734, -1.731, -4.169   , -4.169, 7.667, 2.033],\
      [11.50 , -0.154, -1.557, -3.862   , -3.862, 7.358, 4.491]]

w18 = numpy.array(w_18)
m_13=[[11.22, 1.216 ,  -1.29],\
      [11.00, 1.625 ,  -1.17],\
      [11.00, 1.625 ,  -1.17],\
      [11.00, 1.625 ,  -1.17],\
      [10.87, 1.391 ,  -1.02],\
      [10.87, 1.391 ,  -1.02],\
      [10.81, 1.013 ,  -0.86],\
      [10.81, 0.479 ,  -0.55],\
      [11.03, 0.193 ,  -1.01],\
      [11.49, 0.009 ,  -1.45]]
m13=numpy.array(m_13)      

l13_m   =[10.88,11.03, 11.03, 10.87, 10.71, 10.74, 10.74, 10.74, 10.76, 10.74]
l13_phi =[1.68 , 1.22,  1.22, 2.03 , 1.35 , 0.88 , 0.88 , 0.62 , 0.26 , 0.03 ]
l13_a   =[-0.69,-1.00,-1.00 ,-0.52 ,-0.08 , -0.24, -0.24, -0.22, -0.15, 0.95 ]
l13_phi2=[0.77 , 0.16, 0.16 , 0.29 , 0.67 , 0.33 , 0.33 , 0.15 , 0.14 , 0.09 ]
l13_b   =[-1.42,-1.64,-1.64 ,-1.62 ,-1.46 , -1.6 , -1.6 , -1.6 , -1.6 ,-1.6  ]

d17_m   =[10.78,10.77,10.77 ,10.56 ,10.62 ,10.62 ,10.51 ,10.60 ,10.59 ,10.83 ,11.10 ]
d17_a   =[-1.38,-1.36,-1.36 ,-1.31 ,-1.28 ,-1.28 ,-1.28 ,-1.57 ,-1.67 ,-1.76 ,-1.98 ]
d17_phi =[1.187,1.070,1.070 ,1.428 ,1.069 ,1.069 ,0.969 ,0.295 ,0.228 ,0.090 ,0.016 ]
d17_b   =[-0.43,0.03 ,0.03  ,0.51  ,0.29  ,0.29  ,0.82  ,0.07  ,-0.08 ,0.    ,0.    ]
d17_phi2=[1.92 ,1.68 ,1.68  ,2.19  ,1.21  ,1.21  ,0.64  ,0.45  ,0.21  ,0.    ,0.    ]


#print w18[:,0]
#print w18[:,
fig = figure()
for num in range(1,11):
    #cat=numpy.loadtxt('cos_data/cos_flux_s%s.txt'%(num))
    cat=numpy.loadtxt('cos_data/data_s%s.el'%(num))
    z_min = zbins[int(num) -1]
    z_max = zbins[int(num)]
    print z_min,z_max
    v1=cosmocalc(z_min,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
    v2=cosmocalc(z_max,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
    vmax=(v2-v1)*1.5*sqDeg2sr/(4*pi)

    m=numpy.arange(6.05,13.9,0.1)#0.1    
    m=numpy.arange(6.1,13.8,0.2)#0.2    
    m=numpy.arange(6.25,13.5,0.5)#0.5
    #m=numpy.arange(6.3,13.6,0.6)#0.6
    #m=numpy.arange(6.35,13.5,0.7)#0.7
    #m=numpy.arange(6.4,13.5,0.8)#0.8
    #m=numpy.arange(6.45,13.5,0.9)#0.9
    #m=numpy.arange(6.5,13.5,1)#1
    
    dm=1./(m[1] - m[0])
    
    mf_13=phi(m, l13_m[num-1], l13_phi[num-1]*1e-3, l13_phi2[num-1]*1e-3, l13_a[num-1], l13_b[num-1])
    mf_17=phi(m, d17_m[num-1], d17_phi[num-1]*1e-3, d17_phi2[num-1]*1e-3, d17_a[num-1], d17_b[num-1])
    mf_14=phi(m,m13[:,0][num-1],m13[:,1][num-1]*1e-3, 0.               , m13[:,2][num-1], 0.)
    mf_18=phi(m,w18[:,0][num-1],10**w18[:,3][num-1],10**w18[:,4][num-1], w18[:,1][num-1], w18[:,2][num-1])
    
    fig.add_subplot(4,3,num)
    #fig = figure()
    #if num==1:print m
    
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,0.1))
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,0.2))
    mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,0.5))
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,0.6))
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,0.7))
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,0.8))
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,0.9)) #0.8
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=numpy.arange(6,14,1))
    #mf,mfbins=numpy.histogram(cat[:,-2],bins=m,histtype='step')
    
    #hist(cat[:,-2],bins=numpy.arange(6,14,0.5), weights=dm,histtype='step')
    plot(mfbins[1:],mf*dm,ls='steps')
    #hist(mfbins[:-1],mfbins,weights=mf*dm,histtype='step')
    plot(m,mf_13*vmax,linewidth=3,label='Ilbert+2013')
    plot(m,mf_14*vmax,'-.',linewidth=3,label='Muzzin+2013')
    plot(m,mf_17*vmax,'--',linewidth=3,label='Davidzon+2017')
    #splot(m,mf_18*vmax,'.',markersize=12,label='Wright+2018')
    xlim(6,14)
    yscale('log')
    #print mfbins
    print 'modelb bins  13   17  18'
    co_13 = mf*dm/(vmax*mf_13)
    co_17 = mf*dm/(vmax*mf_17)
    co_18 = mf*dm/(vmax*mf_18)
    
    for i in range(len(m)):
        print '%4.2f %4.2f %6.3f %6.3f %6.3f '%(m[i],mfbins[i],co_13[i]*100,co_17[i]*100,co_18[i]*100)
    
    
    if num==8 or num==9:
    	xticks(numpy.arange(7,14,1),fontsize=22)
    if num == 10:
    	xticks(numpy.arange(7,14,1),fontsize=22)
    if num<4:
        xticks( visible=False)
        #ylim(0,2550)
        ylim(1,5e5)
        text(10,100000,r'$\rm{%.2f < z < %.2f}$ '%(z_min,z_max),alpha=1,fontsize=22 )
        if num==1:
            yticks([1e1,1e2,1e3,1e4,1e5], fontsize=20)
            continue
    elif num<7:
        xticks( visible=False)
        #ylim(0,2800)
        ylim(1,5e5)
        text(10,100000,r'$\rm{%.2f < z < %.2f}$ '%(z_min,z_max),alpha=1,fontsize=22 )
        if num==1 or num ==4:
            yticks([1e1,1e2,1e3,1e4,1e5], fontsize=20)
            continue
    elif num <10:
        ylim(1,5e5)
        text(10,100000,r'$\rm{%.2f < z < %.2f}$ '%(z_min,z_max),alpha=1,fontsize=22 )
        if num==7:
            yticks([1e1,1e2,1e3,1e4,1e5], fontsize=20)
            xticks( visible=False)
            continue
    if num==10:
       #if num==1:
       #ylim(0,440)
       #yticks([0,10,100,1000,10000], fontsize=20)
       text(10,100000,r'$\rm{%.2f < z < %.2f}$ '%(z_min,z_max),alpha=1,fontsize=22 )
       #legend(frameon=False).draggable()
       #else:
       ylim(1,5e5)
       yticks([1e0,1e1,1e2,1e3,1e4,1e5],fontsize=20)
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
    ylim(1,5e5)

legend(frameon=False).draggable()
fig.text(0.06,0.5,r'$\rm{N}$',fontsize=28, ha='center',va='center', rotation='vertical')   
fig.text(0.5,0.06,r'$\rm{Stellar}$ $\rm{Mass (\log_{10}[M/M_\odot])}$',fontsize=28, ha='center',va='center')
#for a
subplots_adjust(hspace=0,wspace=0)
show()#savefig(BOUT_HISTO)

