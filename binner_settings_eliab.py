"""
Parameter file and inits for lumfunc.py and simulate.py
"""

import os,sys,numpy

#-------------------------------------------------------------------------------
# Which program is loading this settings file?
exe=sys.argv[0]
context='cmdline'
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
nb= 40#40#38#40#39#37#41#50#39#37 #13 #24 #27 #34 #37  # 38 or 41
outdir='chains_150414a' # based on 140123a
run_num=outdir.split('_')[-1]
if context=='s' or context=='i': outdir='sims/%s' % outdir.split('_')[-1]
logfile='README.txt'
variablesfile='variables.txt'
comment='Run portability tests'

# What to fit - SPL, TPL, XPL or QPL
nlaws=4

dataset='sdss'
run_num_run='150414a'

# Specify the data file (within outdir) and set up some survey parameters
if dataset=='sdss':
    datafile='sdss_dr12s2.txt'
    SURVEY_AREA=14555 # sq. deg.
    SURVEY_NOISE=150.0 # uJy

#-------------------------------------------------------------------------------

# Parameters for binning a catalogue
if context=='b':
    BIN_CAT_CLIP=None
    CORR_RESOLUTION=None
    CORR_BINS=None
    BOUT_HISTO=os.path.join(dataset,'flux_histo_%s.pdf'%('2'))
    if 'sdss' in dataset:
        BIN_CAT_FORM=6
        BIN_CAT_CLIP=None # ? Cut on something
        CORR_BINS=numpy.ones(nb)
        BIN_COL=3 # - pixel flux
        BIN_CAT=os.path.join(dataset,'dr12_s1.txt')
        CORR_RESOLUTION=1.0
        BOUT_CAT=os.path.join(dataset,'sdss_dr12s1.txt')
        
#-------------------------------------------------------------------------------

# Set up the binninga
binstyle='sdss'

if binstyle=='sdss':
    bins=numpy.linspace(-1000.0,1000.0,41)
    bins=[ -700., -340,  -295,  -251,  -206,\
      -161. ,-138., -116., -93. ,-71.,\
      -48.,   -26. ,  -4. ,18., 40.,   63.26530612, 85.,  108., 130,  153,\
      175,  197., 219 , 242., 264 ,  287., 309,  332., 354.,\
         377.55102041,   422.44897959,   467.34693878,   512.24489796,\
         557.14285714,     646.93877551,\
         736.73469388,   826.53061224,  1006.12244898, \
         1140.81632653 , 1320.40816327,  1500.]        
    bins=numpy.array(bins)
    assert(len(bins)-1==nb)

nbins=len(bins)
dbins=numpy.gradient(bins)

# Warning messages etc.

print 'MOTD: %s' % MOTD

#-------------------------------------------------------------------------------

