"""
Parameter file and inits for lumfunc.py and simulate.py
"""

import os,sys,numpy

#-------------------------------------------------------------------------------
# Which program is loading this settings file?
exe=sys.argv[0]
context='cmdline'
num = sys.argv[1]
if    'reconstruct' in exe: context='r'
elif  'simulate'    in exe: context='s'
elif  'lumfunc'     in exe: context='l'
elif  'plot'        in exe: context='p'
elif  'binner'      in exe: context='b'
elif  'extractor'   in exe: context='e'
elif  'resample'    in exe: context='u'
elif  'inject'      in exe: context='i'
print 'Context is %s' % context
#-------------------------------------------------------------------------------

# Master parameters
MOTD=''
RESUME=False # Turn checkpointing on 
nb=31#21 #22#25#55#40#38#40#39#37#41#50#39#37 #13 #24 #27 #34 #37  # 38 or 41
outdir='chains_150414a' # based on 140123a
run_num=outdir.split('_')[-1]
if context=='s' or context=='i': outdir='sims/%s' % outdir.split('_')[-1]
logfile='README.txt'
variablesfile='variables.txt'
comment='Run portability tests'

dataset='cosmos'
run_num_run='150414a'
dataf = 'cos_s%s.txt'%num
print dataf
#dataf = 'skads_s1_3_noisy.el'
#dataf = 'dr12s_f4.txt'
num =3
# Specify the data file (within outdir) and set up some survey parameters
if dataset=='sdss':
	#if len(dataf)==15: #vol
	print len(dataf)
	#if len(dataf)==15: # 12  16
	#	num = dataf[-9]# vol
		#num = dataf[-5]# -5  -9
	#	print num
	#else:
	#	num = dataf[-6:-4]
	if int(num) > 7:
		SURVEY_AREA=8609.67 # sq. deg.
	else:
		SURVEY_AREA=6672 # dr7 sq. deg.
    
	datafile='sdss_dr12s%s.txt'%(num)
	SURVEY_NOISE=150.0 # uJy
if dataset=='skads': #skads	
	if len(dataf)==17: #skads
		num = dataf[7]# vol
	else:
		num = dataf[7:9]
	
	print num,len(dataf)
	print 'noting happend!'
	#sys.exit()
	SURVEY_AREA=16.
	datafile='data_cos_s%s.txt'%(num)
	SURVEY_NOISE=3.0 # uJy
	#dataset='sdss'

if dataset=='cosmos':
	if len(dataf)==10: 
		num = dataf[5]
	else:
		num = dataf[5:7]
	
	print num,len(dataf)
	print 'noting happend!'
	#sys.exit()
	SURVEY_AREA=1.38
	datafile='data_cos_s%s.txt'%(num)
	SURVEY_NOISE=2.4 # uJy
	#dataset='sdss'

#-------------------------------------------------------------------------------
def get_bins(flux_min, flux_max, z, Lbinwidth=0.4, Lmin=20):
        dl = get_dl(z)
        Lmax = numpy.log10(get_Lbins([flux_max],z,dl)[0])
        Lbins = numpy.arange(Lmin,round(Lmax+1), Lbinwidth)
        positive_sbins = get_sbins(10**Lbins,z,dl)*1e6
        n = numpy.where(positive_sbins > abs(flux_min))[0][0] +1
        #negative_sbins = 0 - positive_sbins[abs(flux_min) >=positive_sbins]
        negative_sbins =numpy.sort( 0 - positive_sbins[:n])
        sbins = numpy.concatenate([negative_sbins,positive_sbins])
        return sbins


# Parameters for binning a catalogue
if context=='b':
    BIN_CAT_CLIP=None
    CORR_RESOLUTION=None
    CORR_BINS=None
    BOUT_HISTO=os.path.join(dataset,'flux_histo_%s.pdf'%(num))
    if 'sdss' in dataset:
        BIN_CAT_FORM=3 #9 #3
        BIN_CAT_CLIP=None # ? Cut on something
        CORR_BINS=numpy.ones(nb)
        BIN_COL=-1 # - pixel flux
        BIN_CAT=os.path.join(dataset,dataf)
        CORR_RESOLUTION=1.0
        BOUT_CAT=os.path.join(dataset,datafile)
        
    elif 'skads' in dataset:
        BIN_CAT_FORM=9 #9 #3
        BIN_CAT_CLIP=None # ? Cut on something
        CORR_BINS=numpy.ones(nb)
        BIN_COL=2 # - pixel flux
        BIN_CAT=os.path.join('cos_data',dataf)
        CORR_RESOLUTION=1.0
        BOUT_CAT=os.path.join('cos_data',datafile)
        
    elif 'cosmos' in dataset:
        BIN_CAT_FORM=9 #9 #3
        BIN_CAT_CLIP=None # ? Cut on something
        CORR_BINS=numpy.ones(nb)
        BIN_COL=-1 # - pixel flux
        BIN_CAT=os.path.join('cos_data',dataf)
        CORR_RESOLUTION=1.0
        BOUT_CAT=os.path.join('cos_data',datafile)
        
#-------------------------------------------------------------------------------

# Set up the binninga
binstyle='sdss'

if binstyle=='sdss':
    #bins=numpy.linspace(-1000.0,1000.0,41)
    #old bin style

    bins =[ -3.65833012e+01,\
        -1.45640745e+01,  -5.79806251e+00,  -2.30825026e+00,\
        -3.65833012e-01,   3.65833012e-01,   9.18930980e-01,\
         2.30825026e+00,   5.79806251e+00,   1.45640745e+01,\
         3.65833012e+01,   9.18930980e+01,   2.30825026e+02,\
         5.79806251e+02,   1.45640745e+03,   3.65833012e+03,\
         9.18930980e+03,   2.30825026e+04,   5.79806251e+04,\
         1.45640745e+05,   3.65833012e+05,   9.18930980e+05]
         
    bins=[-6.73224346e+01,  -2.68015439e+01,\
        -1.06698868e+01,  -4.24775845e+00,  -1.69106310e+00,\
        -6.73224346e-01,  -2.68015439e-01,  -1.06698868e-01,\
        -4.24775845e-02,  -1.69106310e-02,   6.73224346e-03,\
         1.69106310e-02,  4.24775845e-02,   1.06698868e-01,\
         2.68015439e-01,   6.73224346e-01,   1.69106310e+00,\
         4.24775845e+00,   1.06698868e+01,   2.68015439e+01,\
         6.73224346e+01,   1.69106310e+02,   4.24775845e+02,\
         1.06698868e+03,   2.68015439e+03,   6.73224346e+03,\
         1.69106310e+04,   4.24775845e+04,   1.06698868e+05,\
         2.68015439e+05,   6.73224346e+05,   1.69106310e+06]
         
    print 'I want the length'
    print len(bins)-1,nb
    nb = len(bins) -1
    assert(len(bins)-1==nb)
    print bins

nbins=len(bins)
dbins=numpy.gradient(bins)

# Warning messages etc.

print 'MOTD: %s' % MOTD

#-------------------------------------------------------------------------------

