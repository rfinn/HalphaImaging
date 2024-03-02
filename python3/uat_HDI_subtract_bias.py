#!/usr/bin/env python

'''
GOAL


PROCEDURE


REQUIRED MODULES
ccdproc


'''

import glob
import os
import sys
from astropy.io import fits
import argparse
import astropy

from astropy.nddata import CCDData
import ccdproc
from astropy.modeling import models

parser = argparse.ArgumentParser(description ='Subtract overscan bias and trim images.')
parser.add_argument('--filestring', dest = 'filestring', default = 'tc', help = 'string to use to get input files (default = "c", which grabs all of the files "c*o00.fits")')
parser.add_argument('--verbose', dest='verbose', default=False,action='store_true', help='print extra messages for troubleshooting')

#parser.add_argument('--gain', dest = 'gain', default= 1.3, help = 'gain in e-/ADU.  default is 1.3, which applies to HDI camera')
#parser.add_argument('--rdnoise', dest = 'rdnoise', default= 7.3, help = 'gain in e-/ADU.  default is 1.3, which applies to HDI camera')
args = parser.parse_args()

# get list of bias frames
biasfiles = glob.glob(args.filestring+'*b00.fits')

# combine bias
if args.verbose:
    print(f"bias images = {biasfiles}")
if len(biasfiles) < 3:
    print('Not enough bias frames')
    sys.exit()

############################################################        
# combine bias images using average combine, sigma clip
############################################################

combinedbias = ccdproc.combine(biasfiles,method='average',sigma_clip=True,unit=u.adu,\
                       sigma_clip_low_thresh=5, 
                       sigma_clip_high_thresh=5,
                       sigma_clip_func=np.ma.median,
                       sigma_clip_dev_func=astropy.stats.mad_std)

combinedbias.data = combinedbias.data.astype('f4')
combinedbias.header['NBIAS'] = nbias
combinedbias.header['COMBINED'] = True

combinedbiasfile = 'bias_average.fits'
print('  Writing {}'.format(combinedbiasfile))
combinedbias.write(combinedbiasfile, overwrite=True)

############################################################
# mv bias frames to BIAS subdirectory
############################################################
if not os.path.exists('BIAS_FRAMES'):
    os.mkdir('BIAS_FRAMES')
for b in biasfiles:
    os.rename(b,'BIAS_FRAMES/'+b)


############################################################
# subtract bias from other images
############################################################

# get list of bias frames
flatfiles = glob.glob(args.filestring+'*f00.fits')
# get list of bias frames
sciencefiles = glob.glob(args.filestring+'*o00.fits')

allfiles = flatfiles + sciencfiles


for f in allfiles:
    # subtract the bias, and pre-pend a 'b'

    # read in image
    print('working on ',f)
    # convert data to CCDData format and save header
    ccd = CCDData.read(f, unit='adu')

    # subtract overscan
    head_updates = {'SUBBIAS':True}
    b_subtracted = ccdproc.subtract_bias(ccd, bias, add_keyword = head_updates)
    
    b_subtracted.write('b'+f, overwrite=True)    
    
