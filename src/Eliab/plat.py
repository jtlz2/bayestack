import numpy
import os,sys,math,shutil
from pylab import*
from matplotlib.ticker import AutoMinorLocator
from cosmocalc import cosmocalc
from scipy.interpolate import interp1d
z_c, c = numpy.loadtxt('completeness', unpack=1)
C_z= interp1d(z_c,c)

Ho=70
wm=0.3
fmin = 1e-29
specindx = 0.7

Zall,E = numpy.loadtxt('../bayestack/z_E.txt',unpack=True)
def openlit(name):
    f=open(name)
    lines=f.readlines()
    zslice=[]
    data=[]
    if 'pine' in name:
        f2=open('mcalpine-2013-sfg.txt')
        lines2=f2.readlines()
    print 'start'
    for i in range(len(lines)):
        line=lines[i]
        lin=line.split()
        if line[0]=='#':
           if line[1]=='z':
            data.append(numpy.array(zslice))
            zslice=[]
           
           #print line
           continue
           
        #print line
        if 'pine' in name:
            #line2=lines2[i].split()
            err= numpy.sqrt(float(lin[2])**2)# + float(line2[2])**2)
            phi=float(lin[1])#+float(line2[1])
            log_err = 0.434*err/phi
            zslice.append([numpy.log10(float(lin[0])), numpy.log10(phi) , log_err])
        else:
            zslice.append([float(lin[0]), float(lin[1]), float(lin[2]), float(lin[3])])
    data.append(numpy.array(zslice))
    return numpy.array(data)

def open_anytype(name,cols,zr=[],F='mJy',L=True,S=False,getz=False,Mi=False,ug=False,band='L'):
	if Mi:
		cols.append(Mi)
		if ug:
		  cols.append(ug)
		  cols.append(ug+1)
		  z,s,mi,u,g = numpy.loadtxt(name,usecols=cols,unpack=True)
		else:
		  z,s,mi = numpy.loadtxt(name,usecols=cols,unpack=True)
	else:
		z,s = numpy.loadtxt(name,usecols=cols,unpack=True)
	if zr==[]:
		pass
	else:
		s = s[numpy.where((z> zr[0]) & (z< zr[1]) )]
		if Mi:
			mi=mi[numpy.where((z> zr[0]) & (z< zr[1]) )]
		z = z[numpy.where((z> zr[0]) & (z< zr[1]) )]
	
	
	if S:
		if Mi:
	         return s[s>0],z[s>0],dl[s>0],mi[s>0]
		return s 
	if F=='Jy':
	    s*= 1e-26
	elif F=='mJy':
	    s*= 1e-29
	else: #uJy
	    s*= 1e-32
	    
	dl = [cosmocalc(zi,H0=Ho,WM = wm)['DCMR_Mpc'] for zi in z]
	dl =numpy.array(dl)* 3.08e22
	l = math.pi*4 *s* dl**2 * (1 + z)**(specindx+1)
	if band=='S':
	    l*=1#(1.4/3)**(-.7)#before 27Feb this was just l*=1!
	if L:
	    if ug:return l[l>0],z[l>0],mi[l>0],u[l>0],g[l>0]
	    if Mi:return l[l>0],z[s>0],mi[s>0]
	    if getz: return l[l>0],z[s>0]
	    return l[l>0]
	print z[0],s[0],l[0],dl[0]
	print len(l)
	return calcu_zeff(z,s,dl,l)	

#def calcu_zeff(z,s,d,l):
def calcu_zeff(z,l,fmin= 1e-32):
	#E  = dl**2 * (1+z)**(0.7-1)
	l = numpy.array(l)
	ef = l/(4*math.pi*fmin)*3.08
	zef = numpy.interp(ef,E,Zall)
	return zef


def calc_Vmax(zef,l,bins,zup,zlo,skyarea,getLF=True,table=False):
	Vmax =[]
	dz=0.1

	if not type(zef)==float:
	    z  = [min(zup,ze) for ze in zef]
	else:
	    z = 0.5*(zup+zlo)*(1+ numpy.zeros(len(l)))
	    print 'here?'
	
	V1 = [cosmocalc(zi,H0=Ho,WM = wm)['VCM_Gpc3']*1e9 for zi in z]
	V2 = cosmocalc(zlo,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
	#zmean = (z + zlo)/2
	V1 = numpy.array(V1)
	vmax_old = skyarea*(V1-V2)/(4*pi)
	
	
	
	for i in range(len(z)):
	 zup=z[i]
	 z_i=zlo 
	 vmax =0 
	 while z_i < zup:
	    
	    V1 = cosmocalc(z_i,H0=Ho,WM = wm)['VCM_Gpc3']*1e9 
	    z_i+=dz
	    V2 = cosmocalc(z_i,H0=Ho,WM = wm)['VCM_Gpc3']*1e9
	
	    vmax += skyarea*(V2-V1)/(4*pi)
	 Vmax.append(vmax)
	 #print vmax_old[i],vmax
	
	#for i in range(len(vmax)):
	 #if vmax[i]<0:
	 # print z[i],zef[i],zup
	  #sys.exit()
	#for i in range(len(l)):
	#print 'done with Vmax old ',numpy.mean(vmax_old),' vs new ',numpy.mean(Vmax)
	Vmax=numpy.array(Vmax)
	#sys.exit()
		
		
	#	Vmax.append(abs(vmax))
	if getLF:
	   return lumfunc(l,bins,Vmax,table)
	else:
	   return Vmax#vmax_old#Vmax

def lumfuncs(N,Lbins,Vmax,return_dm=False):
    dm = [(Lbins[i+1] - Lbins[i])/0.4 for i in range(len(Lbins)-1)]
    dm.append(dm[-1])
    #print dm [-3]
    #sys.exit()
    if return_dm:
      return N/Vmax/dm, dm[-8]
    return N/Vmax/dm

def lumfunc(L,bins,Vmax,table=False):
	bin_norm = numpy.zeros(len(bins))
	err_bin  = numpy.zeros(len(bins))
	rho_m 	 = numpy.zeros(len(bins))
	bin_size = numpy.zeros(len(bins))
	
	L = log10(numpy.array(L))
	o_Vmax    =  1./numpy.array(abs(Vmax))
	o_logVmax = log10(o_Vmax)

	index = numpy.digitize(L,bins = bins,right=0)
	
	#print len(index),len(Vmax), len(err_bin), len(o_Vmax)
	#print L
	#print index
	#sys.exit()
	#bins, his = numpy.histogram(L,bins = bins)#,weights = o_Vmax)
	for i in range(len(Vmax)):
		indx = index[i]
		#print indx
		if indx>=len(bins):
			#print '\n\n\n\n\n'
			#print 'Index error! check boss'
			#print '\n\n\n\n\n'
			continue
		err_bin[indx] += (o_Vmax[i])**2
		bin_norm[indx]+= o_Vmax[i]
		bin_size[indx]+=1
		#print indx
		#bins[indx]	  += 1.
		
	#dm = 0.4/(bins[1] - bins[0]) #binwidth, this should be 0.4
	dm = [(bins[i+1] - bins[i])/0.4 for i in range(len(bins)-1)]
	dm.append(dm[-1])
	#print bin_norm
	#sys.exit()
	print 'Start'
	bin_mid=[]

	for i in range(len(bins)-1):
		try:
		    rho_m[i] = log10(bin_norm[i]/(dm[i]) )
		except:
		    rho_m[i]=0
		print round((bins[i]+ bins[i+1])/2.,2),'    ', int(bin_size[i])
		bin_mid.append((bins[i]+ bins[i+1])/2.)
		#print round(bins[i],2), round(bins[i+1],2), round((bins[i]+ bins[i+1])/2.,2), rho_m[i], bin_size[i]


	#bins = [log10(weight)	for (weight) in bin_norm]
	err_bin = err_bin**0.5
	lin_rhom = 10**rho_m #back to linear space
	#error = [(er*bin1 )/10**(bin1)  for (err,bin1) in zip(err_bin,bins)]
	error = [abs(log10(rho + err)) - abs(log10(rho)) for(err,rho) in zip(err_bin,lin_rhom)]
		
	log_err = 0.434*err_bin/lin_rhom
	
	if table:
	    return rho_m,log_err, bin_mid,bin_size
	return rho_m,log_err
'''
#bins = arange(23,28,0.8)
bins = arange(20.8,24.,0.1)
#bins = bins = arange(24.6,29.2,0.4)
zlo = 1.8
zup = 2.5
area1 =  2.66 #NVSS
area2 =  (3212* 3.0462e-4) #FIRST DR9
area3 = (11663* 3.0462e-4)#FIRST DR7
area4 = (8000* 3.0462e-4)#FIRST DR12
#zef1,l1 = open_multiS()
#rho_1,er1 = calc_Vmax(zef1,l1,bins,zup,zlo,area4)
L = open_flux()
rho,err = calc_Vmax(99.,L,bins,0.45,0.2,1.)

fig,ax = subplots()
errorbar(bins,rho,yerr=err,fmt = 'o')
#errorbar(bins,rho_1,yerr=er1,fmt = 'o',color= 'k',label='NVSS-SDSS DR7 Condon et al 2013')
xlabel('1.4 GHz log[L(W/Hz)]',fontsize = 18)
ylabel(r'log[$\rho_m$(mpc$^{-3}$mag$^{-1}$)]',fontsize = 18)
tick_params(axis='both',which = 'major', labelsize=15,width =2)
tick_params(axis='both',which = 'minor', labelsize=12, width=1)
ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())

#show()


zef1,l1 = opencondon('condons1.txt')
zef2,l2 = open_first('dr12_1d_18.el') 
zef3,l3 = open_first2('matchedr12_1.txt')
rho_1,er1 = calc_Vmax(zef1,l1,bins,zup,zlo,area1)
rho_2,er2 = calc_Vmax(zef2,l2,bins,zup,zlo,area4)
rho_3,er3 = calc_Vmax(zef3,l3,bins,zup,zlo,area4)

#subplot(121)
fig,ax = subplots()
errorbar(bins,rho_2,yerr=er2,fmt = 'o',label='FIRST-SDSS DR12')
errorbar(bins,rho_1,yerr=er1,fmt = 'v',color= 'r',label='Condon et al 2013')
#errorbar(bins,rho_3,yerr=er3,fmt = 'x',color= 'r',label='FIRST')
xlabel('1.4 GHz log[L(W/Hz)]',fontsize = 18)
ylabel(r'log[$\rho_m$(mpc$^{-3}$mag$^{-1}$)]',fontsize = 18)
xlim(25,29.8)
ylim(-12.8,-8)
text(26.5,-10.5,'1.8 < z < 2.5 \n m$_r$ < 18.5',fontsize =16,alpha = 0.6)
tick_params(axis='both',which = 'major', labelsize=15,width =2)
tick_params(axis='both',which = 'minor', labelsize=12, width=1)

ax.xaxis.set_minor_locator(AutoMinorLocator())
ax.yaxis.set_minor_locator(AutoMinorLocator())
legend()
show() 	

bins = arange(22.4,27.6,0.4)
zlo = 0.2
zup = 0.45
zef1,l1 = opencondon('condons2.txt')
zef2,l2 = open_first('matchede3s2.txt')
zef3,l3 = open_first2('matchedr7fir_2.txt')
rho_1,er1 = calc_Vmax(zef1,l1,bins,zup,zlo,area1)
rho_2,er2 = calc_Vmax(zef2,l2,bins,zup,zlo,area2)
rho_3,er3 = calc_Vmax(zef3,l3,bins,zup,zlo,area3)

#subplot(122)
errorbar(bins,rho_2,yerr=er2,fmt = 'o',label='FIRST-SDSS DR9')
errorbar(bins,rho_1,yerr=er1,fmt = 'v',color= 'r',label='NVSS-SDSS DR7 Condon et al 2013')
errorbar(bins,rho_3,yerr=er3,fmt = 'x',color= 'k',label='FIRST-SDSS DR7 ')
xlabel('1.4 GHz log[L(W/Hz)]',fontsize = 18)
ylabel(r'log[$\rho_m$(mpc$^{-3}$mag$^{-1}$)]',fontsize = 18)
xlim(18,28)
ylim(-10.8,-3.8)
text(23.6,-9.4,'0.2 < z < 0.45 \n $M_i$ < -23',fontsize = 14 )
tick_params(axis='both',which = 'major', labelsize=15,width =2)
tick_params(axis='both',which = 'minor', labelsize=12, width=1)
legend()
show() 
'''
