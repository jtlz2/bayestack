from pylab import*
from astropy.io import fits
from astropy import wcs
import numpy,sys
from utils import gaussian
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit


data = fits.getdata('/home/eliab/Documents/OC/Project/cosmos/vla_3ghz_msmf.fits')
h = fits.getheader('/home/eliab/Documents/OC/Project/cosmos/vla_3ghz_msmf.fits')
w = wcs.WCS(h)

bmaj=h['BMAJ']
bmaj=0.75
A=numpy.pi*bmaj*bmaj/(4*numpy.log(2))
pix_l=0.2
flux_d_to_flux=A/(pix_l**2)
factor=2.6 # for 1x1 grid (3x3)
print 'Beam area ',A
print 'pixel area ',pix_l**2
print 'flux_d to flux',flux_d_to_flux
z=[0.1,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,4.0]
#z=[0.1,0.3,0.4,0.6,0.8,1.0,1.3,1.6,2.0,2.5,3.2,3.6,4.0]
num=int(sys.argv[-1])
print z[num-1], ' < z < ', z[num]
        

#ra_d,dec_d,flux_d,npix,multi = numpy.loadtxt('/home/eliab/Documents/OC/Project/Cosmos/VLA_detected/table1.dat',usecols=(3,5,7,-3,-1),unpack=1)
'''
def get_pix_size(real_size_fwhm,z):
    ''
    This function calculates the pixel size of a sources
    given it's FWHM axis major axis and redshit
    it uses the sizes in SFG_size
    ''
    real_size=/2.3548*2*numpy.sqrt(-numpy.log(0.01)) # i get sigma, then calc the diameter where the height is 1percent of the total height
    dl_dpc=numpy.array(get_dL(z))*1e3 #using kpc because the size in in kpc
    ang_s = ((1+z)**2*ang_size)/(dl_dpc)*206264.81 #using angular dist with dL
    return numpy.sqrt(ang_s**2+0.75**2)/0.2 #convloving size with beam
'''
def myg(x,a,b,norm):
    return norm* gaussian(x,a,b,norm=False)

def get_total_flux(data,c_x,c_y,r_x=3,r_y=3):
    flx=0.
    c_x=int(c_x)
    c_y=int(c_y)
    
    for i in range(c_x-r_x,c_x+1+r_x):
        for k in range(c_y-r_y,c_y+1+r_y):
            flx+= data[0][0][i][k]
            
    return flx/flux_d_to_flux
                      
'''
r_x,r_y = 1,1
print 'extracting from ', r_x*2+1, 'x',r_y*2+1, ' box'
flux_x =numpy.zeros(len(ra_d))
for i in range(len(flux_x)):
    r,d,o,p=w.wcs_world2pix(ra_d[i],dec_d[i],0,1.501191666680E+02,2.205833333330)
    #flux_x[i] = data[0][0][float(d)-1][float(r)-1]
    flux_x[i] = get_total_flux(data,float(d)-1,float(r)-1,r_x,r_y)

#f = open('/home/eliab/Documents/OC/Project/Cosmos/cos_d.txt','w')

#f.write()
#f.write('%10s %10s %12s %12s'%('Ra','Dec','flux','xflux\n'))
#for i in range(len(flux_x)):
#   f.write('%10f %10f %12f %12f \n'%(ra_d[i],dec_d[i],flux_d[i],flux_x[i]*1e6))
#f.close()
#print 'done with main cos'
#sys.exit()

x=logspace(1,5,100)

plot(flux_d,2.6*flux_x*1e6,'o',fillstyle='none')
plot(flux_d[npix>100],2.6*flux_x[npix>100]*1e6,'sm', label='Npixes > 100')
plot(flux_d[multi==1],2.6*flux_x[multi==1]*1e6,'og', label='Multi-component')
plot(x,x,'-r',linewidth=3)
xlabel('Intergated flux $\mu$Jy',fontsize=25)
ylabel(' Extrated flux $\mu$Jy',fontsize=25)

xscale('log')
yscale('log')
xlim(7,1e5)
ylim(7,1e5)
tick_params(axis='both',which = 'minor', labelsize=10, width=2)
tick_params(axis='both',which = 'major', labelsize=20,width =3)
legend(frameon=False,numpoints=1, prop={'size': 20}).draggable()
show()
'''
cos_lim = numpy.loadtxt('cos_data/data_s%s.el'%num)
ra,dec = cos_lim[:,0],cos_lim[:,1]
ID = cos_lim[:,-1]
Ra =ra + ((-0.041 * ra + 6.059)/3600.)
Dec=dec + ((0.058 *dec - 0.147)/3600.)
Ra=ra
Dec=dec
print len(cos_lim)

surv, IDs, flx =[],[],[]

f=open('cosmos_counter.el')
lines = f.readlines()
for line in lines:
    line=line.split()
    if '#' in line:continue
    surv.append(line[3])
    IDs.append(float(line[4]))
    flx.append(float(line[8]))
surv=numpy.array(surv)
IDs=numpy.array(IDs)
flx=numpy.array(flx)    
    
flux = numpy.zeros(len(cos_lim))
for i in range(len(flux)):
    r,d,o,p=w.wcs_world2pix(Ra[i],Dec[i],0,1.501191666680E+02,2.205833333330)
    flux[i] = get_total_flux(data,float(d)-1,float(r)-1)
    
fluxes =numpy.array(flux)*1e6
print 'done extracting'
f = open('cos_data/cos_s%s.txt'%num,'w')
f2 = open('cos_data/cos_flux_s%s.txt'%num,'w')
c=1
#f.write()
myflux=[]
taflux=[]
for i in range(len(fluxes)):
   flux=fluxes[i]
   if ID[i] in IDs:
        if surv[IDs==ID[i]] =='COSMOS2015':
            myflux.append(fluxes[i])
            if flx[IDs==ID[i]][0]>500:
                
                if flux>50:
                    flux = flx[IDs==ID[i]][0]
                    
                print fluxes[i], flx[IDs==ID[i]][0], 'using ',flux
            f2.write('%15f %15f %7f %10f %10f \n'%(Ra[i],Dec[i],cos_lim[i][2],fluxes[i],flx[IDs==ID[i]][0]))
            taflux.append(flx[IDs==ID[i]][0])
            #print 'it works', i+1, c
            c=c+1
   f.write('%15f %15f %7f %10f \n'%(Ra[i],Dec[i],cos_lim[i][2],flux))
f.close()

SURVEY_NOISE = 2.5
binwidth=0.4*SURVEY_NOISE
#fig,ax = plt.subplots()
plot(taflux,myflux,'o')

ylabel(r'$\rm{Extracted}$ $\rm{flux density (\mu Jy)}$', fontsize=25)
xlabel(r'$\rm{Catalogue}$ $\rm{flux density (\mu Jy)}$',fontsize=25)
#xlabel(r'$\rm{Smolcic}$ $\rm{et}$ $\rm{al}$ $\rm{2017}$$\rm{flux-densities}$',fontsize=25)
plot(taflux,taflux,'r' )
xlim(4,1e4)
ylim(4,1e4)
xscale('log')
yscale('log')
text(200,8, '%.2f < z < %.2f'%(z[num-1],z[num]),fontsize=20)
#ax.yaxis.set_minor_locator(AutoMinorLocator())
#ax.xaxis.set_minor_locator(AutoMinorLocator())
tick_params(axis='both',which = 'minor', labelsize=10, width=2)
tick_params(axis='both',which = 'major', labelsize=20,width =3)

fig,ax = plt.subplots()

n,b,p = hist(fluxes, bins=numpy.arange(min(fluxes),(20.0*SURVEY_NOISE)+binwidth,binwidth),histtype='step',color='blue')
norm=max(n)
coef, unc = curve_fit(myg,b[1:],n,p0=[0.1,1,500])
print coef, unc
y =norm* gaussian(b,0.6,SURVEY_NOISE*1.45,norm=False)
y2 = numpy.max(n)*gaussian(b,0,SURVEY_NOISE*1.45,norm=False)
y3 =coef[2]* gaussian(b,coef[0]-5.29755559e-01,coef[1],norm=False)
#plot(b,y2,'r--',linewidth=2,label=r'$\rm{Guassian(\mu=0)}$')
plot(b,y3,'--r',linewidth=2,label=r'$\rm{Guassian}$ $\rm{fit}$')
ylabel(' Number of sources',fontsize=25)
xlabel(' Extrated flux $\mu$Jy',fontsize=25)
axvline(coef[1]
*5,color='g',linewidth=3)
xlim(-12,17)
ax.yaxis.set_minor_locator(AutoMinorLocator())
ax.xaxis.set_minor_locator(AutoMinorLocator())
tick_params(axis='both',which = 'minor', labelsize=10, width=2)
tick_params(axis='both',which = 'major', labelsize=20,width =3)
text(-5,20, '%.2f < z < %.2f'%(z[num-1],z[num]),fontsize=20)
#legend(frameon=False,numpoints=1, prop={'size': 20}).draggable()
show()

'''
f_x,f_d,sep = numpy.loadtxt('/home/eliab/Documents/OC/Project/Cosmos/cosmos_d%s_0_3arc.txt'%num,unpack=True, usecols=[7,6,-1])
x=logspace(1,5,100)

plot(f_d,f_x,'o',fillstyle='none')
#plot(flux_d[npix>100],flux_x[npix>100]*1e6,'sm', label='Npixes > 100')
#plot(flux_d[multi==1],flux_x[multi==1]*1e6,'og', label='Multi-component')
plot(x,x,'-r',linewidth=3)
xlabel('Intergated flux $\mu$Jy',fontsize=25)
ylabel(' Extrated flux $\mu$Jy',fontsize=25)

xscale('log')
yscale('log')
xlim(7,1e3)
ylim(7,1e3)
tick_params(axis='both',which = 'minor', labelsize=10, width=2)
tick_params(axis='both',which = 'major', labelsize=20,width =3)
legend(frameon=False,numpoints=1, prop={'size': 20})#.draggable()



for i in range(1,6):
    f_x,f_d,sep = numpy.loadtxt('/home/eliab/Documents/OC/Project/Cosmos/cosmos_d%s_0_3_%sarc.txt'%(num,i),unpack=True, usecols=[7,6,-1])

    fig,ax = plt.subplots()
    plot(f_d,f_x,'o',fillstyle='none')
#plot(flux_d[npix>100],flux_x[npix>100]*1e6,'sm', label='Npixes > 100')
#plot(flux_d[multi==1],flux_x[multi==1]*1e6,'og', label='Multi-component')
    plot(x,x,'-r',linewidth=3)
    xlabel('Intergated flux $\mu$Jy',fontsize=25)
    ylabel(' Extrated flux $\mu$Jy',fontsize=25)

    xscale('log')
    yscale('log')
    xlim(7,1e3)
    ylim(7,1e3)
    tick_params(axis='both',which = 'minor', labelsize=10, width=2)
    tick_params(axis='both',which = 'major', labelsize=20,width =3)
    legend(frameon=False,numpoints=1, prop={'size': 20})#.draggable()

 
show()
''
bash
xt to get cosmos_d"$i"_0_3_"$j"arc.txt ; ./stilts tskymatch2 ifmt1=ascii ifmt2=ascii in1=Cosmos/cos_s$i.txt in2=Cosmos/cos_d_$j.txt ofmt=asci out=Cosmos/cosmos_d"$i"_0_3_"$j"arc.txt ra1=col1 ra2=Ra dec1=col2 dec2=Dec error=0.3 join=1and2 find=best ; done; done

for i in 3 4 5 6 7 8 9 10; do echo using Cosmos/cos_s$i.txt and Cosmos/cos_d.txt to get cosmos_d"$i"_0_3arc.txt ; ./stilts tskymatch2 ifmt1=ascii ifmt2=ascii in1=Cosmos/cos_s$i.txt in2=Cosmos/cos_d.txt ofmt=asci out=Cosmos/cosmos_d"$i"_0_3arc.txt ra1=col1 ra2=Ra dec1=col2 dec2=Dec error=0.3 join=1and2 find=best ; done
'''
