"""
Parameter file and inits for bayestack.py, lumfunc.py and simulate.py
"""

import os,sys,numpy,math
from utils import sqDeg2sr,beamFac
from cosmocalc import cosmocalc

#-------------------------------------------------------------------------------
# Which program is loading this settings file?
exe=sys.argv[0]
context='cmdline'
if    'reconstruct' in exe: context='r'
elif  'simulate'    in exe: context='s'
elif  'lumfunc'     in exe: context='l'
elif  'bayestack'   in exe: context='y'
elif  'plot'        in exe: context='p'
elif  'binner'      in exe: context='b'
elif  'extractor'   in exe: context='e'
elif  'resample'    in exe: context='u'
elif  'inject'      in exe: context='i'
elif  'dnds_lumfunc'in exe: context='z'
print 'Context is %s' % context
#-------------------------------------------------------------------------------

# New-style settings <-- bayestack.py

#dataset='video'
binStyle=1
nlaws=3 # Number of ppl laws (or poly coefficients) to fit
floatNoise=True
modelFamily='LFlognorm_dpl'#'ppl' 'poly' # 'LFsch' # 'LFdpl'
outdir='chains_201102i' # based on 140123a
noiseZones=[0] # Select which noise zones to analyse
doPoln=False
doRayleigh=False # Set Rice signal to zero, i.e. Rayleigh distribution
doRedshiftSlices=True # For joint z slice analysis
whichRedshiftSlices=[9]
setBreaks=False
Breaks=[450]
# LF variables
LF_SPEC_INDEX=0.7
Ho=70.0
wm=0.3

#-------------------------------------------------------------------------------

# Master parameters
#MOTD=''
RESUME=True # Turn checkpointing on
#nb= 22#40#38#40#39#37#41#50#39#37 #13 #24 #27 #34 #37  # 38 or 41
dnds0=False # Leave as False otherwise no recon line...
binsHigh=False # Run with this True then set to False
#outdir='chains_150514a' # based on 140123a
run_num=outdir.split('_')[-1]
#if context=='l': outdir=os.path.join('/home/jtlz2/bartolomeu/output',outdir)
if context=='s' or context=='i': outdir='sims/%s' % outdir.split('_')[-1]
logfile='README.txt'
variablesfile='variables.txt'
comment='SDSS QSO RLF - tests'
SEED_SAMP=9123 # [-1 for clock]
breaks = [] #breaks in flux
#floatNoise=False

dataset='cosmos'


# Specify the data file (within outdir) and set up some survey parameters
numRedshiftSlices=len(whichRedshiftSlices)
numNoiseZones=5

if dataset == 'first':
    datafile='ketron2013_jz.txt'
    SURVEY_AREA=2.16 # sq. deg.
    SURVEY_NOISE=150.0 # uJy
elif dataset == 'mca':
    datafile='bondi2003_mca.txt'
    SURVEY_AREA=1.00
    SURVEY_NOISE=17.0
elif dataset=='sdss':
    SURVEY_NOISE=3.0 # uJy
    #SURVEY_AREA=8609.67 # sq. deg.
    SURVEY_AREA=8. #Dr7
    zmanifestf='sdss/sdss_manifest.txt'
    zmanifest=numpy.genfromtxt(zmanifestf)
    datafiles=['sdss/sdss_dr12s%i.txt'%i for i in zmanifest[[s-1 for s in whichRedshiftSlices],0]]
    if len(datafiles[0])==20: # 12  16
                num = datafiles[0][-5]# -5  -9
    else:
                num = datafiles[0][-6:-4]

    if int(num) > 7:
		SURVEY_AREA=8. # sq. deg.
    else:
		SURVEY_AREA=8. # dr7 sq. deg.
    #datafiles=['sims/161101d/sim.txt']
    #datafile= datafiles[0]
    #numRedshiftSlices=zmanifest.shape[0]
    redshifts=zmanifest[[s-1 for s in whichRedshiftSlices],3]
    numRedshiftSlices=len(whichRedshiftSlices)

elif dataset=='cosmos':
    SURVEY_NOISE=2.3 # uJy
    SURVEY_AREA=1.45
    zmanifestf='cos_data/cos_manifest.txt'
    zmanifest=numpy.genfromtxt(zmanifestf)
    datafiles=['cos_data/data_cos_s%i.txt'%i for i in zmanifest[[s-1 for s in whichRedshiftSlices],0]]

    if len(datafiles[0])==10: 
		num = datafiles[0][5]
    else:
		num = datafiles[0][5:7]

    redshifts=zmanifest[[s-1 for s in whichRedshiftSlices],3]
    numRedshiftSlices=len(whichRedshiftSlices)


#-------------------------------------------------------------------------------

# Command-line arguments [mostly vestigial]
verbose=False
sigma=SURVEY_NOISE # uJy
area=SURVEY_AREA # sq. deg.
#model=1 # Vary all parameters
#loud=False
#tellthetruth=False


#-------------------------------------------------------------------------------

###REF_AREA=10.0 # sq. deg. for following arrays - vestigial - ***keep this fixed

# LF parameters

# Prior types
LMIN_PRIOR='U'
LMAX_PRIOR='U'
LMIN2_PRIOR='U'
LMAX2_PRIOR='U'
LSTAR_PRIOR='U'
LSLOPE_PRIOR='U'
LNORM_PRIOR='U'
LZEVOL_PRIOR='U'
LSLOPE2_PRIOR='U'
LSLOPE2_2_PRIOR='U'

#LOGNORM
LSTAR_2_PRIOR='U'
LSLOPE_2_PRIOR='U'
LNORM_2_PRIOR='U'
LSIGMA_PRIOR='GAUSS'

# Prior ranges
LMIN_MIN=18
LMIN_MAX=23.3	
LMAX_MIN=27.
LMAX_MAX=30

LMIN2_MIN=19.2
LMIN2_MAX=21.6 
LMAX2_MIN=26
LMAX2_MAX=30

LSTAR_MIN=23.5
LSTAR_MAX=27.5
LSLOPE_MIN=-0.5
LSLOPE_MAX=2
LNORM_MIN=-9
LNORM_MAX=-4.7
LSLOPE2_MIN=-5
LSLOPE2_MAX=1


LSTAR_2_MIN=20.0
LSTAR_2_MAX=24.9
LSLOPE_2_MIN= -5.0
LSLOPE_2_MAX=5.0
LSLOPE2_2_MIN= -5.0
LSLOPE2_2_MAX=5.0

LNORM_2_MIN=-5.5
LNORM_2_MAX=2
LSIGMA_MIN=0.6
LSIGMA_MAX=0.1

if setBreaks:
    dl = cosmocalc(redshifts[0],H0=Ho,WM=wm)['DL_Mpc']
    l=math.pi*4 *(Breaks[0]*1e-32)* (dl* 3.08e22)**2 * (1 + redshifts[0])**((LF_SPEC_INDEX)+1)
    LSTAR_MIN=LSTAR_MAX=23.1
    #LMAX_MIN=LMAX_MAX=24.6
    print 'setting break'
    LSTAR_PRIOR = 'DELTA'
    #LMAX_PRIOR  = 'DELTA'
    
    

# LF redshift evolution (if > 1 z slice only)
LZEVOL_MIN=-5.0
LZEVOL_MAX=5.0
if len(whichRedshiftSlices) <= 1:
    LZEVOL_MIN=LZEVOL_MAX=0.0 # delta fn on LZEVOL

#S0_MIN=S0_MAX=S0_TRUE            # delta fn on S0
#BETA_MIN=BETA_MAX=BETA_TRUE    # delta fn on beta
#D_MIN=D_MAX=D_TRUE            # delta fn on C

NOISE_MIN=0.5*SURVEY_NOISE
#NOISE_MAX=2.0*SURVEY_NOISE
NOISE_MAX=2.0*SURVEY_NOISE
if not floatNoise:
    NOISE_MIN=NOISE_MAX=SURVEY_NOISE # delta fn on NOISE

 #		triangle plot settings
PLOT_LABEL=''
triangle='triangle_%s.png' % run_num

#-------------------------------------------------------------------------------

# Experiment resolution (for source-confusion prior/s)
confusionNoiseThreshold=0.5 # x sigma
checkConfusionNoise=False
confusionNoiseCheckSmax=5.0*SURVEY_NOISE # or None => Smax for each itern

checkSmin=False
numericalCheckSmin=True
dOmega=7.72e-6 # This is 10 x 10 arcsec^2 / sq. deg.
#dOmega=7.72e-7 # This is 10 arcsec^2 / sq. deg.
#dOmega=2.44e-7 
#dOmega=2.235e-10
#dOmega=2.7778e-3
# Confusion-noise factor
div=2.0
# Non-overlapping prior on Smin/max

# Set some MultiNEST parameters here
outstem='1-'

#n_live_points=1000 # was 1000
multimodal=True
max_modes=3
SEED_SAMP=SEED_SAMP # [-1 for clock]


# Switch to INS
do_INS=False
n_live_points=500

max_iter=0
evidence_tolerance=0.5 # was 0.5
sampling_efficiency=0.3 # default is 0.8
