#!/usr/bin/env python

"""
GOAL
* run on each image from the NOAO pipeline
* subtract sky from each ccd
* run getzp on each ccd
* calculate offset for each amp
* scale the image and ivar accordingly
* save sky-subtracted, scaled image

USAGE:

python BOK_pipline_fixampoffsets.py image_name image_filter

"""

import os
import sys

homedir = os.getenv("HOME")
sys.path.append(homedir+'/github/HalphaImaging/python3/')


from astropy.io import fits
import numpy as np

from subtract_median import subtract_median
import getzp

image_name = sys.argv[1]
ivar_name = image_name.replace('ooi','oow')
image_filter = sys.argv[2]


### AMPS

NAXIS1 = 4032
NAXIS2 = 4096


class args():
    """ for replicating argparse input to getzp """
    def __init__(self,image,instrument,filter,nexptime=False):
        self.image = image
        self.instrument = instrument
        self.filter = filter
        self.normbyexptime = nexptime
        self.verbose = False
        self.nsigma = 3.5
        
        self.d = os.getenv("HOME")+'/github/HalphaImaging/astromatic/'
        self.catalog = None
        self.fwhm = None
        self.useri = False
        self.mag = 0
        self.naper = 5
        self.fit = False
        self.flatten = 0
        self.norder = 2
        

#################################################################
### MAIN PROGRAM
#################################################################

# subtract sky from each ccd
files = [image_name]
if not os.path.exists('m'+image_name):
    subtract_median(files,overwrite=False,MEF=True)
else:
    print(f"Found median-subtracted image for {image_name} - using this.")

# run getzp on each ccd
hdu = fits.open('m'+image_name) # read in median-subtracted image
ihdu = fits.open(ivar_name) # read in median-subtracted image
for h in range(1,len(hdu)):
    hdu[h].header.set('EXPTIME',hdu[0].header['EXPTIME'])
    hdu[h].writeto(f'temp{h}.fits',overwrite=True)
    medsubimage = 'm'+image_name
    myargs = args(f'temp{h}.fits','i',image_filter,nexptime=True)    
    zp = getzp.getzp(myargs)

    zp.getzp()


    # calculate offset for each amp
    #    residual_all = 10.**((magfit - yplot)/2.5)        
    #    self.residual_all = residual_all
    #    plt.scatter(self.matchedarray1['X_IMAGE'][self.fitflag],self.matchedarray1['Y_IMAGE'][self.fitflag],c = (residual_all),vmin=v1,vmax=v2,s=15)
    overall_med = np.ma.median(zp.residual_all)
    print(f"Overall median of residuals = {overall_med:.3f}")
    # calculate the offset for each amplifier and scale the data accordingly
    for ix in range(2):
        xmin = 0 + NAXIS1//2*ix
        xmax = NAXIS1//2 * (ix +1)

        for iy in range(2):
            ymin = 0 + NAXIS2//2*iy
            ymax = NAXIS2//2 * (iy+1)
            print(xmin,xmax,ymin,ymax)
            flag = (zp.residual_allx > xmin) & (zp.residual_allx < xmax) &\
                (zp.residual_ally > ymin) & (zp.residual_ally < ymax) 
            amp_med = np.ma.median(zp.residual_all[flag])
            print(f"median for quadrant {ix+ iy} = {amp_med:.3f}")
            amp_scale = overall_med/amp_med

            # scale the data
            # scale the image and ivar accordingly            
            hdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*hdu[h].data[ymin:ymax,xmin:xmax]

            # what is the correct way to scale the weights?  same as image or inverse???
            # in the weight image, high numbers are good
            ihdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*ihdu[h].data[ymin:ymax,xmin:xmax]
            


# save sky-subtracted, scaled image
hdu.writeto('zm'+image_name,overwrite=True)

ihdu.writeto('zm'+ivar_name,overwrite=True)
