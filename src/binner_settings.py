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
    datafile='sdss_dr12s1.txt'
    SURVEY_AREA=400.0 # sq. deg.
    SURVEY_NOISE=140.0 # uJy

#-------------------------------------------------------------------------------

# Parameters for binning a catalogue
if context=='b':
    BIN_CAT_CLIP=None
    CORR_RESOLUTION=None
    CORR_BINS=None
    BOUT_HISTO=os.path.join(dataset,'flux_histo_%s.pdf'%(run_num))
    if 'sdss' in dataset:
        BIN_CAT_FORM=6
        BIN_CAT_CLIP=None # ? Cut on something
        CORR_BINS=numpy.ones(nb)
        BIN_COL=3 # - pixel flux
        BIN_CAT=os.path.join(dataset,'dr12_s1_nonzero.txt')
        CORR_RESOLUTION=1.0
        BOUT_CAT=os.path.join(dataset,'sdss_dr12s1.txt')
        
#-------------------------------------------------------------------------------

# Set up the binning
binstyle='sdss'

if binstyle=='sdss':
    #bins=numpy.linspace(-270.0,1000.0,8)
    bins=numpy.array([-80.0,-50.0,-20.0,-10.0,-5.0,0.0,5.0,10.0,20.0,\
                        50.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0])
    bins=numpy.linspace(-1000.0,1000.0,41)
    #bins=numpy.array([-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,85.0])
    #bins[0]=-69.0
    assert(len(bins)-1==nb)

nbins=len(bins)
dbins=numpy.gradient(bins)

# Warning messages etc.

print 'MOTD: %s' % MOTD

#-------------------------------------------------------------------------------

