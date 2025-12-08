#!/usr/bin/env python


'''

GOAL: subtract median from images

OVERVIEW:
* designed to complete final stages of coaddition for INT WFC images
* 


USAGE:
* to run on a list of image

python ~/github/HalphaImaging/python3/subtract_sky.py --mef --filestring1 ksb --filestring2 r_v1.fits

* to run on one image

python ~/github/HalphaImaging/python3/subtract_sky.py --mef –oneimage ksb_220428_040758_ooi_r_v1.fits

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
import math
print(homedir)

sys.path.append('/home/rfinn/github/HalphaImaging/python3/')
#import imutils
import ha_imutils as imutils
import ccdproc as ccdp
import glob
from astropy.io import fits
from astropy import stats
import numpy as np


##########################################################
### SUBTRACT MEDIAN FROM IMAGE
##########################################################

def subtract_median(files,overwrite=False,MEF=False,MOS=False):
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
        if 'weight' in fname:
            continue
        if 'bpm' in fname:
            continue
        subtract_median_one(fname,MEF=MEF,MOS=MOS)

        
def subtract_median_one(fname,MEF=True,overwrite=False,MOS=False):
    if not overwrite:
        if os.path.exists("m"+fname) and not overwrite:
            print("m"+fname,' already exists.  moving to next file')
            return
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
            d,median,std = imutils.subtract_median_sky(hdu[i].data.copy(),getstd=True)

            #print('median for hdu {} = {}'.format(i,median))
            #print('check if median is nan: {}'.format(median == np.nan))
            #print('check if median == nan: {}'.format(median == nan))
            if (str(median) == 'nan'):                    
                print('using alternate median for hdu ',i)
                mmean, mmed,mstd = stats.sigma_clipped_stats(hdu[i].data,sigma=3)
                print('alternate estimate of median = {:.2f}'.format(mmed))
                if mmed is not np.nan:
                    hdu[i].data -= mmed
                    median = mmed
                    hdu[i].header.set('MEDSUB',value=mmed,comment='sky med subtract_med')                 
                    hdu[i].header.set('SKYMED',value=mmed,comment='sky med subtract_med')
                    hdu[i].header.set('SKYSTD',value=mstd,comment='sky std subtract_med')
                
            else:
                hdu[i].data = d
                hdu[i].header.set('MEDSUB',value=median,comment='sky med subtract_med')                
                hdu[i].header.set('SKYMED',value=median,comment='sky med subtract_med')
                hdu[i].header.set('SKYSTD',value=std,comment='sky std subtract_med')                
            
    elif MOS:
        # get weight image
        weight_imname = fname.replace('.coadd','.coadd.weight')
        whdu = fits.open(weight_imname)
        
        ymax, xmax = hdu[0].data.shape 
        # case for mosaic data, subtract median in each CCD/AMP
        # these are boundaries that Becky measured from ds9
        xvals = [1, 1018, 2068, 3138, 4190, 5265, 6330, 7383, xmax]
        yvals = [1, 4130, ymax]

        xvals = np.array(xvals,'i')
        yvals = np.array(yvals,'i')
        # subtract 1 to convert to python zero indexed system
        xvals = xvals - 1
        yvals = yvals - 1
        
        buffer = 10 # buffer for measuring sky

        # keep track of overall statistics
        
        average_med = 0
        average_std = 0
        namp = 1
        for j in range(len(xvals)-1):
            xmin = xvals[j]
            xmax = xvals[j+1]

            for k in range(len(yvals) -1 ):
                ymin = yvals[k]
                ymax = yvals[k+1]
            
                d,median,std = imutils.subtract_median_sky(hdu[0].data[ymin+buffer:ymax-buffer,xmin+buffer:xmax-buffer],getstd=True,subtract=False,\
                                                               weightimage = whdu[0].data[ymin+buffer:ymax-buffer,xmin+buffer:xmax-buffer])

                average_med += median
                average_std += std
                # subtract median
                hdu[0].data[ymin:ymax,xmin:xmax] -= median
                hdu[0].header.set('REGION'+str(namp),value=f"{xmin}:{xmax},{ymin}:{ymax}",comment='{xmin}:{xmax},{ymin}:{ymax}')                               
                hdu[0].header.set('MEDSUB'+str(namp),value=median,comment='sky med subtract_med')                
                hdu[0].header.set('SKYMED'+str(namp),value=median,comment='sky med subtract_med')
                hdu[0].header.set('SKYSTD'+str(namp),value=std,comment='sky std subtract_med')
                namp += 1
        # calculate average med and std
        average_med = average_med/(namp-1)
        average_std = average_std/(namp-1)
        hdu[0].header.set('SKYMED',value=average_med,comment='ave sky med subtract_med all amps')
        hdu[0].header.set('SKYSTD',value=average_std,comment='ave sky std subtract_med all amps')
        
    else:
        # background subtraction
        hdu[0].data,median,std = imutils.subtract_median_sky(hdu[0].data,getstd=True)
        if math.isnan(median):
            print(f"WARNING: could not subtract median for {fname}")
        else:
            hdu[0].header.set('MEDSUB',value=median,comment='median subtraction')
            hdu[0].header.set('SKYMED',value=median,comment='sky med subtract_med')
            hdu[0].header.set('SKYSTD',value=std,comment='sky std subtract_med')

    
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
    parser.add_argument('--mos',  default = False, action='store_true', help = 'set this for MOS files')    
    parser.add_argument('--oneimage', dest = 'oneimage', default = None,help = 'supply an image name to run sky subtraction on one image')    
    # TODO
    # add a mosaic flag, where median will be calculated for each amplifier (8 ccds, and each ccd has two amplifiers)
    # determine sky in each amp and subtract
    # need pixel box of full amplifier and pixel box to use for calcu
    args = parser.parse_args()

    #if args.hdi:
    #    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']
    #else:
    #    keys = ['naxis1', 'naxis2', 'imagetyp', 'filter', 'exptime','instrmnt']
    if args.oneimage is not None:
        subtract_median_one(args.oneimage,overwrite=args.overwrite,MEF=args.mef, MOS=args.mos)
    else:
        matchstring = args.filestring+'*.fits'
        if args.filestring2 is not None:
            matchstring = args.filestring+'*'+args.filestring2
        files = glob.glob(matchstring)
        files.sort()
        #print(files)
        subtract_median(files,overwrite=args.overwrite,MEF=args.mef, MOS=args.mos)
