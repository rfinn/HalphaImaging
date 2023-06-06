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

python BOK_pipline_fixampoffsets.py image_name 

"""

import os
import sys
import shutil

homedir = os.getenv("HOME")
sys.path.append(homedir+'/github/HalphaImaging/python3/')


from astropy.io import fits
import numpy as np

from subtract_median import subtract_median
import getzp

image_name = sys.argv[1]
ivar_name = image_name.replace('ooi','oow')

if image_name.find('r_v1') > -1:
    image_filter = 'r'
if image_name.find('Ha+4nm') > -1:
    image_filter = 'ha'
if image_name.find('Ha4nm') > -1:
    image_filter = 'ha'
dq_name = image_name.replace('ooi','ood')

### AMPS

NAXIS1 = 4032
NAXIS2 = 4096


class args():
    """ for replicating argparse input to getzp without using argparse"""
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
        self.nofixbok = True

#################################################################
### MAIN PROGRAM
#################################################################

# subtract sky from each ccd
files = [image_name]
if not os.path.exists('m'+image_name):
    subtract_median(files,overwrite=False,MEF=True)
else:
    print()
    print(f"Found median-subtracted image for {image_name} - using this.")
    print()
    
# run getzp on each ccd
hdu = fits.open('m'+image_name) # read in median-subtracted image
ihdu = fits.open(ivar_name) # read in median-subtracted image
image_name_base = image_name.replace('.fits','')
allresid = []
allresid1d =  []
allresidx = []
allresidy = []
allzp = []
ccd_medians = np.zeros(len(hdu)-1)
firstpass = True
print("")
#print(f"RA = {hdu[0].header['CRVAL1']:.6f}
for h in range(1,len(hdu)):
    ##
    # this is stepping through each CCD and running getzp separately on it
    # it then stores the ZP of each image
    # and will scale the images so that each CCD has the same ZP
    #
    # but this doesn't account for the variations in each amplifier, does it?
    # that is taken care of in the next section.
    # here we just save the x and y positions and residuals so that amp
    # adjustment can be made below.
    ##
    hdu[h].header.set('EXPTIME',hdu[0].header['EXPTIME'])
    hdu[h].header.set('OBJECT',hdu[0].header['OBJECT'])
    hdu[h].header.set('FILTER',hdu[0].header['FILTER'])
    ##
    # write out the individual ccd image to use as input to getzp
    ##
    hdu[h].writeto(f'{image_name_base}-temp{h}.fits',overwrite=True)
    medsubimage = 'm'+image_name

    ##
    # set up arguments to pass into getzp
    # using instrument=i b/c we don't have the full image yet
    # for getzp to help normalize
    ##

    # NOTE:
    # changing instrument from 'i' to 'b' so that color corrections are appropriate
    # I added a flag in getzp to allow us to skip the function that adjust for amp offsets
    # in the coadded image, as opposed to the MEF images that we are working with here
    myargs = args(f'{image_name_base}-temp{h}.fits','b',image_filter,nexptime=True)
    print(f'running getzp on ccd {h}')
    zp = getzp.getzp(myargs)

    zp.getzp()
    allresid.append(zp.residual_all)
    allresid1d.extend(zp.residual_all)    
    allresidx.append(zp.residual_allx)
    allresidy.append(zp.residual_ally)
    allzp.append(-1*zp.zp)
##
# get global median for all ccds, so that ccds are normed relative to each other
##
zp_ref = np.mean(np.array(allzp))
global_med = np.median(allresid1d)
print(f"global median for all ccds = {global_med:.3f}")
for h in range(1,len(hdu)):
    ##
    # calculate the offset for each amplifier and scale the data accordingly
    ##
    quad = 0
    print(f"CCD {h}:")
    ccd_med = np.ma.median(allresid[h-1])
    zp_scale = 10.**((zp_ref - allzp[h-1])/2.5)
    print(f"zp ccd = {allzp[h-1]:.2f}, zpref = {zp_ref:.2f}, zpscale = {zp_scale:.3f}")
    for ix in range(2):
        xmin = 0 + NAXIS1//2*ix
        xmax = NAXIS1//2 * (ix +1)

        for iy in range(2):
            ymin = 0 + NAXIS2//2*iy
            ymax = NAXIS2//2 * (iy+1)
            #print(xmin,xmax,ymin,ymax)
            flag = (allresidx[h-1] > xmin) & (allresidx[h-1] < xmax) &\
                (allresidy[h-1] > ymin) & (allresidy[h-1] < ymax) 
            amp_med = np.ma.median(allresid[h-1][flag])

            amp_scale = ccd_med/amp_med * zp_scale

            print(f"\tmedian and scale for quadrant {quad} = {amp_med:.3f} {amp_scale:.3f}")
            ##
            # scale the data so that each ccd/amp has the same ZP
            # scale the image and ivar accordingly
            ##
            hdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*hdu[h].data[ymin:ymax,xmin:xmax]

            ##
            # what is the correct way to scale the weights?  same as image or inverse???
            # in the weight image, high numbers are good
            ##
            ihdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*ihdu[h].data[ymin:ymax,xmin:xmax]
            quad += 1


# save sky-subtracted, scaled image
hdu.writeto('zm'+image_name,overwrite=True)

ihdu.writeto('zm'+ivar_name,overwrite=True)


shutil.copy(dq_name,'zm'+dq_name)


# check to see if header files exist
# if they do, then copy to have z prefix

scamp_header = 'm'+image_name.replace('.fits','.head')
if os.path.exists(scamp_header):
    shutil.copy(scamp_header,'z'+scamp_header)


# clean up temp files
#for i in range(1,5):
#    os.remove(f'{image_name_base}-temp{i}.fits')
#    os.remove(f'n{image_name_base}-temp{i}.fits')    
#    os.remove(f'n{image_name_base}-temp{i}.cat')
    #os.remove(f'n{image_name_base}-temp{i}_pan_tab.csv')    
