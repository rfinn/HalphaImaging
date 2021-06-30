#!/usr/bin/env python


'''

GOAL: subtract median from images

OVERVIEW:
* designed to complete final stages of coaddition for INT WFC images
* 

REFERENCES:
https://photutils.readthedocs.io/en/stable/background.html

https://reproject.readthedocs.io/en/stable/mosaicking.html

ccdproc wcs projection
https://ccdproc.readthedocs.io/en/latest/image_combination.html
'''

#from photutils import make_source_mask
import os
homedir = os.getenv("HOME")
import sys
print(homedir)

sys.path.append('/home/rfinn/github/HalphaImaging/python3/')
import imutils
import ccdproc as ccdp
import glob
from astropy.io import fits
from astropy import stats
import numpy as np


##########################################################
### SUBTRACT MEDIAN FROM IMAGE
##########################################################

def subtract_median(files,overwrite=False,MEF=False):
    '''
    INPUT:
    * files - list of files for median subtraction

    OPTIONAL INTPUT:
    * overwrite - overwrite the original image with the median subtracted image; default is False; when false, a new images is created with m prepended to the filename
    * MEF - flag to indicate images are multi-extension format, like for 90prime; default is False.

    OUTPUT:
    * the function creates median-subtracted images of each file in filelist
    * median-subtracted images have "m" pre-pended to input image name

    '''
    print('subtracting median from images')
    for fname in files:
        if not overwrite:
            if os.path.exists("m"+fname) and not overwrite:
                print("m"+fname,' already exists.  moving to next file')
                continue
            else:
                print("{} -> m{}".format(fname,fname))
        else:
            print("{} -> {}".format(fname,fname))

        # read in image 
        hdu = fits.open(fname,memmap=False)
 
        if MEF:
            # if MEF flag is set, assume primary header is extension 0
            # loop over additional extenstions and subtract median
            nextensions = len(hdu)
            for i in range(1,nextensions):
                d,median = imutils.subtract_median_sky(hdu[i].data.copy())
                hdu[i].data = d
                print('median for hdu {} = {}'.format(i,median))
                #print('check if median is nan: {}'.format(median == np.nan))
                #print('check if median == nan: {}'.format(median == nan))
                if (str(median) == 'nan'):                    
                    print('using alternate median for hdu ',i)
                    mmean, mmed,mstd = stats.sigma_clipped_stats(hdu[i].data,sigma=3)
                    print('alternate estimate of median = {:.2f}'.format(mmed))
                    if mmed is not np.nan:
                        hdu[i].data -= mmed
                        median = mmed
                        hdu[i].header.set('MEDSUB',value=median,comment='median subtraction')
                
                else:
                    hdu[i].header.set('MEDSUB',value=median,comment='median subtraction')
            
            
        else:
            # background subtraction
            hdu[0].data,median = imutils.subtract_median_sky(hdu[0].data)
            if median is not np.nan:
                hdu[0].header.set('MEDSUB',value=median,comment='median subtraction')
        if overwrite:
            hdu.writeto(fname,overwrite=True)
        else:
            hdu.writeto("m"+fname,overwrite=True)
        hdu.close()
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description ='Subtract the median from images, after masking out objects and growing mask.')

    parser.add_argument('--filestring', dest = 'filestring', default = 'WFC', help = 'filestring to match. default is WFC')
    parser.add_argument('--filestring2', dest = 'filestring2', default = None, help = 'second filestring to match. default is None.  set to ooi for 90prime data.')    
    parser.add_argument('--overwrite', action = 'store_true', default = False, help = 'overwrite file?  the default is false, so that a new file with m prefix is created.')
    parser.add_argument('--mef', action = 'store_true', default = False, help = 'set this for MEF files, like with 90prime')        
    args = parser.parse_args()

    #if args.hdi:
    #    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']
    #else:
    #    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']

    matchstring = args.filestring+'*.fits'
    if args.filestring2 is not None:
        matchstring = args.filestring+'*'+args.filestring2
    files = glob.glob(matchstring)
    files.sort()
    #print(files)
    subtract_median(files,overwrite=args.overwrite,MEF=args.mef)
