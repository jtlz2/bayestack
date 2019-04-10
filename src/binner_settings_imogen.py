"""
Parameter file and inits for binner.py
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
nb= 13#40#38#40#39#37#41#50#39#37 #13 #24 #27 #34 #37  # 38 or 41
outdir='chains_150526b' # based on 140123a
run_num=outdir.split('_')[-1]
if context=='s' or context=='i': outdir='sims/%s' % outdir.split('_')[-1]
logfile='README.txt'
variablesfile='variables.txt'
comment='Run portability tests'

# What to fit - SPL, TPL, XPL or QPL
nlaws=4

#dataset='10C_LH'
dataset='10C_LH_t2'
run_num_run='150526b'

# Specify the data file (within outdir) and set up some survey parameters
if dataset=='sdss':
    datafile='sdss_dr12s1.txt'
    SURVEY_AREA=400.0 # sq. deg.
    SURVEY_NOISE=140.0 # uJy
elif dataset == '10C_LH':
    datafile='10C_LH_binned.txt' #I think this should be binned sc file, not flux list
    SURVEY_AREA=0.156 # sq. deg. 
    SURVEY_NOISE=21.0 # uJy
elif dataset == '10C_LH_t2':
    datafile='10C_LH_binned_t2.txt' 
    SURVEY_AREA=0.310 # sq. deg. 
    SURVEY_NOISE=45.0 # uJy

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
    elif dataset == '10C_LH':
        BIN_CAT_FORM=5
        BIN_CAT='10C_LH/WSRT_pix_vals.txt' #Lockman hole deep source catalogue
        BIN_CAT_CLIP=None
        BIN_COL=0 #best flux = col 15 for 10C cat in uJy
        BOUT_CAT='10C_LH/10C_LH_binned_wider.txt'
    elif dataset == '10C_LH_t2':
        BIN_CAT_FORM=7
        BIN_CAT='10C_LH_t2/WSRT_pix_vals_t2.txt' #Lockman hole deep source catalogue
        BIN_CAT_CLIP=None
        BIN_COL=0 #best flux = col 15 for 10C cat in uJy
        BOUT_CAT='10C_LH_t2/10C_LH_binned_t2.txt'
        
#-------------------------------------------------------------------------------

# Set up the binning
#binstyle='10C_LH'
binstyle='10C_LH_t2'

if binstyle=='sdss':
    #bins=numpy.linspace(-270.0,1000.0,8)
    bins=numpy.array([-80.0,-50.0,-20.0,-10.0,-5.0,0.0,5.0,10.0,20.0,\
                        50.0,100.0,200.0,300.0,400.0,500.0,600.0,700.0,800.0])
    bins=numpy.linspace(-1000.0,1000.0,41)
    #bins=numpy.array([-80.0,-65.0,-50.0,-40.0,-30.0,-20.0,-15.0,-10.0,-8.0,-6.5,-5.0,-3.5,-2.5,-1.0,-0.65,-0.5,-0.25,-0.1,-0.05,0.0,0.05,0.1,0.25,0.5,0.65,1.0,2.5,3.5,5.0,6.5,8.0,10.0,15.0,20.0,30.0,40.0,50.0,65.0,85.0])
    #bins[0]=-69.0
    assert(len(bins)-1==nb)

elif binstyle=='10C_LH':
    #original bins
    #bins=numpy.array([-67.0,-60.0,-40.0,-30.0,-20.0,-10.0,-5.0,0.0,5.0,10.0,15.0,20.0,25.0,30.0,35.0,40.0,50.0,60.0,70.0,80.0,110.0])
    #wider bins
    bins=numpy.array([-67.0,-60.0,-20.0,0.0,20.0,40.0,80.0,120.0,200.0])
    bins[0]=-67.0

elif binstyle=='10C_LH_t2':
    bins=numpy.array([-355.0,-140.0,-80.0,-40.0,-20.0,-10.0,0.0,10.0,20.0,40.0,60.0,80.0,120.0,225.0])
    bins[0]=-355.0

nbins=len(bins)
dbins=numpy.gradient(bins)

# Warning messages etc.

print 'MOTD: %s' % MOTD

#-------------------------------------------------------------------------------

