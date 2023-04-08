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
allresid = []
allresid1d =  []
allresidx = []
allresidy = []
allzp = []
ccd_medians = np.zeros(len(hdu)-1)
firstpass = True
for h in range(1,len(hdu)):
    hdu[h].header.set('EXPTIME',hdu[0].header['EXPTIME'])
    hdu[h].writeto(f'temp{h}.fits',overwrite=True)
    medsubimage = 'm'+image_name
    myargs = args(f'temp{h}.fits','i',image_filter,nexptime=True)
    print(f'running getzp on ccd {h}')
    zp = getzp.getzp(myargs)

    zp.getzp()
    allresid.append(zp.residual_all)
    allresid1d.extend(zp.residual_all)    
    allresidx.append(zp.residual_allx)
    allresidy.append(zp.residual_ally)
    allzp.append(-1*zp.zp)
# get global median for all ccds, so that ccds are normed relative to each other
zp_ref = np.mean(np.array(allzp))
global_med = np.median(allresid1d)
print(f"global median for all ccds = {global_med:.3f}")
for h in range(1,len(hdu)):
    # calculate offset for each amp
    #    residual_all = 10.**((magfit - yplot)/2.5)        
    #    self.residual_all = residual_all
    #    plt.scatter(self.matchedarray1['X_IMAGE'][self.fitflag],self.matchedarray1['Y_IMAGE'][self.fitflag],c = (residual_all),vmin=v1,vmax=v2,s=15)


    # calculate the offset for each amplifier and scale the data accordingly
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
            # scale the data
            # scale the image and ivar accordingly            
            hdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*hdu[h].data[ymin:ymax,xmin:xmax]

            # what is the correct way to scale the weights?  same as image or inverse???
            # in the weight image, high numbers are good
            ihdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*ihdu[h].data[ymin:ymax,xmin:xmax]
            quad += 1


# save sky-subtracted, scaled image
hdu.writeto('zm'+image_name,overwrite=True)

ihdu.writeto('zm'+ivar_name,overwrite=True)


shutil.copy(dq_name,'zm'+dq_name)
