#!/usr/bin/env python

"""
Simple wrapper for stacker.py in the context of lumfunc.py
"""

import sys
import importlib
import stacker
from profile_support import profile


__name_cached=__name__
param_file=sys.argv[-1]

#setf='%s' % param_file.split('.')[-2]
# Import the settings variables
if __name__=='__main__':
    setf='settings'
else:
    setf='%s'%param_file

set_module=importlib.import_module(setf)
globals().update(set_module.__dict__)
__name__=__name_cached

#-------------------------------------------------------------------------------

@profile
def main():

    """
    """

    # Read the radio data
    print 'Reading radio map %s' % perGalRadioMap
    imdata,imhead=stacker.readRadio(fradio=perGalRadioMap)
    # Read the noise map
    print 'Reading noise map %s' % perGalNoiseMap
    noisedata,noisehead=stacker.readRadio(fradio=perGalNoiseMap)

    # Read the K-band catalogue (up to z=5)
    # Remember to set the cuts inside stacker.readCatalogue (GENERALIZE THIS)
    gals=stacker.readCatalogue(do_shuffle=False,version=5,filename=perGalMasterCat,\
                               cutKFlags=perGalCutKFlags,chibestthresh=\
                               perGalChiBestThresh,cutHaloFlags=perGalCutHaloFlags,\
                               zReadThresh=perGalzReadThresh,
                               reliabilityReadThresh=perGalReliabilityThresh,\
                               IDRange=perGalIDRange)

    # Set up the D_L curve (do once)
    cosmo=stacker.calculateLumdist()

    # Optionally scramble the redshifts
    if perGalzScramble:
        gals=stacker.scrambleRedshifts(gals,dz=perGalzScrambleDZ,\
                                       seed=perGalzScrambleDZSeed)

    # Use the D_L curve to generate the Dunne cosmology factors
    gals=stacker.updateCosmology(gals,cosmo)

    # Optionally scramble the masses
    if perGalMScramble:
        gals=stacker.scrambleMasses(gals,dm=perGalMScrambleDM,\
                                    seed=perGalMScrambleDMSeed)

    # Optionally shuffle the alphas and deltas independently
    #if perGalAlphaDeltaShuffle:
    #    gals=stacker.shufflePosns(gals)

    # *** Should set the clips and thresholds to infinity in the below?
    if not perGalXYScramble:
        stacker.wrapMzCut(gals,imdata,imhead,noisedata,\
                          noise_thresh=perGalNoiseThresh,\
                          clip_thresh=perGalClipThresh,kind=perGalGalaxyType,\
                          wrap_root=perGalWrapRoot,zbins=perGalzBins,\
                          Mbins=perGalMBins,masses=perGalMasses,\
                          KLIM=perGalKLim,kabsbins=perGalkabsBins,\
                          upsampling_factor=perGalUpsampleRate,\
                          upsampled_map=upsampledMap,\
                          upsampled_noise=upsampledNoise,\
                          use_radio_posn=perGalRadioPosn,\
                          do_sim_posns=do_sim_posns,\
                          scramble_posns=perGalXYRandomize,\
                          noise_only=noise_only,\
                          scramble_posns_seed=perGalXYRandomizeSeed,\
                          scramble_posns_mask=perGalXYRandomizeMask)
    # random - Scramble the positions 10 times
    # *** Should set the mass and z ranges here
    elif perGalXYScramble:
        stacker.realise(gals,imdata,imhead,noisedata,n=perGalXYScrambleN,\
                        m=perGalXYScrambleM,noise_thresh=perGalNoiseThresh,\
                        clip_thresh=perGalClipThresh,root=perGalWrapRoot)

    return 0

#-------------------------------------------------------------------------------

if __name__ == '__main__':

    ret=main()
    sys.exit(ret)


    
