#!/usr/bin/env python

'''
GOAL:
  The goal of this program is to convolve the R and Halpha images to a common FWHM before creating mosaics.


PROCEDURE:
  - create a list of all of the images for a specific target
  - run PSFEx to estimate FWHM of each image
  - read in output from PSFEx
  - look for outliers and remove any images that have bad seeing compared to others
  - get max FWHM for remaining images
  - compute gaussian FWHM needed to bring each image to the FWHM_max
  - convolve images


EXAMPLE:



INPUT/OUPUT:

REQUIRED MODULES:

REQUIRED SOFTWARE:
- PSFEx

EXTRA NOTES:

WRITTEN BY:
Rose Finn, started who knows when

EDITED BY:
Rose Finn, 2024-07-25

'''


#!/usr/bin/env python
import os
import numpy as np 
import glob
from matplotlib import pyplot as plt


from astropy.io import fits
from astropy.table import Table
from astropy.convolution import convolve, convolve_fft, Gaussian2DKernel
from astropy.modeling.models import Gaussian2D

import argparse



def get_fwhm(input_images): #measure FWHM of SE catalogs
    nfiles = len(input_images)
    image_fwhm = np.zeros(nfiles,'f')
    image_fwhm_std = np.zeros(nfiles,'f')   
    for i in range(len(input_images)):
        t = input_images[i].split('.fits')
        se_cat = t[0]+'.cat'
        data = fits.getdata(se_cat,2)
        # select unsaturated stars using : class_star > 0.9 and (10 < m < 13)
        flag = (data.CLASS_STAR > 0.9) & (data.MAG_AUTO > 10.) & (data.MAG_AUTO < 13.)
        image_fwhm[i] = np.mean(data.FWHM_IMAGE[flag])
        image_fwhm_std[i] = np.std(data.FWHM_IMAGE[flag])
    return image_fwhm, image_fwhm_std

def get_fwhm_psfex(psfout=None):
    """
    read in output from psfex

    RETURN:
    images = list of image names
    fwhm = l
    """
    # open file
    #from astropy.io.votable import parse
    #votable = parse(psfout)
    if psfout == None:
        psfout = 'psfex.xml'
    pt = Table.read(psfout,table_id='PSF_Fields')
    images = []
    for im in pt['Catalog_Name']:
        images.append(im.replace('.cat','.fits'))
    fwhm = np.array(pt['FWHM_Mean'])
    return images, fwhm



                    
#All of the above is new  stuff


if __name__ == '__main__':
    #The goal of this program is to have a python-based convolution routine
    parser = argparse.ArgumentParser(description ='This code will convolve image cutouts that have bad focus so that we can get a more precise continuum subtraction.')
    parser.add_argument('--prefix', dest = 'prefix', default = 'pointing-1',  help = 'Input the string of images to be convolved (before continuum subtraction, after cutouts). Enter prefix pointing (e.g. pointing-1)')
    parser.add_argument('--gscale',dest = 'gscale', default = .95,help='scale the max fwhm to get the desired output.  default = 0.95')
    parser.add_argument('--convolve',dest = 'convolve', default = False, action ='store_true',help='set this to convolve images.  by default, psfex is run but we do not run convolution.')
    parser.add_argument('--cthreshold',dest = 'cthreshold', default = 0.5, help='threshold for applying convolution to image.  Images with fwhm within cthreshold pixels of fwhm_max will not have convolution applied. default value is 0.5 pixels. the input image is copied over as g*.')
    parser.add_argument('--filestring', dest = 'filestring', default = 'mh', help = 'string to use to get input files (default = "mh", which grabs all of the files "mh*o00.fits")')    
    args = parser.parse_args()



    #search_prefix = args.string+'*-Ha.fits'
    #search_prefix = args.prefix+'*coadd.fits'
    #input_images = glob.glob(search_prefix)
    if not os.path.exists('PSFEX_OUTPUT'):
        os.mkdir("PSFEX_OUTPUT")
    os.system("cp ~/github/HalphaImaging/astromatic/default.psfex .")

    # run psfex on the images
    os.system('ls '+args.filestring+'*o00.cat > input_list')
    os.system('psfex @input_list -XML_NAME psfex.xml')
    
    
    image_names, image_fwhm = get_fwhm_psfex()
    #print(image_fwhm)
    #Need worst FWHM to convolve to
    fwhm_max=np.max(image_fwhm)

    if not args.convolve:
        # print psfex output for each image

        # or write to file

        pass

        print('#######################################')
        print('COMPARING ORIGINAL FWHM')
        print('#######################################')
        print("image   fwhm    ")
        print("------------------------")
        for i in range(len(image_names)):
            print(f"{image_names[i]}: {image_fwhm[i]:.2f} ")

        # write to a file
        newtab = Table([image_names,image_fwhm])
        newtab.write("psfex_original_fwhm.csv",format='csv',overwrite=True)

    
    
    if args.convolve:
        # keep track of image with max FWHM so we don't run convolution on this
        allindex = np.arange(len(image_names))
        max_fwhm_index = allindex[image_fwhm == fwhm_max][0]
        
        print('the largest FWHM = ',fwhm_max)

        # skip any images that are within of fwhm_max
        convolve_flag = np.abs(image_fwhm - fwhm_max) > float(args.cthreshold)
        
        # convolve all images to worst seeing
        # use pyraf.iraf.gauss (sigma = FWHM/2.35)
        # (sigma_out)^2 = (sigma_in)^2 + (sigma_filter)^2
        #
        # (sigma_filter) = np.sqrt[(fwhm_out/2.35)^2 - (sigma_in/2.35)^2] 
        sigma_filter = np.sqrt((float(args.gscale)*fwhm_max/2.35)**2 - (image_fwhm/2.35)**2)
        #convolve_images()

        for i in range(len(image_names)):
            print("#########################################################")
            print(f"working on image {image_names[i]} with kernel={sigma_filter[i]:.3f}")
            print("#########################################################")
            print()
            #if image_fwhm[i] == fwhm_max:
            #    continue
            hdu = fits.open(image_names[i])[0]

            if not convolve_flag[i]:
                print("\t not running convolution on image with max fwhm ",image_names[i])
                chdu = fits.PrimaryHDU(hdu.data, header=hdu.header)
                chdu.writeto('g'+image_names[i], overwrite=True) 
            else:
                kernel = Gaussian2DKernel(sigma_filter[i])
                convolved_data = convolve(hdu.data, kernel)

                chdu = fits.PrimaryHDU(convolved_data, header=hdu.header)
                chdu.writeto('g'+image_names[i], overwrite=True) 

        # run source extractor on convolved images
        # copy psfex files as default
        os.system('cp ~/github/HalphaImaging/astromatic/default.param.psfex default.param')
        os.system('cp ~/github/HalphaImaging/astromatic/default.sex.HDI.psfex default.sex.HDI')        
        os.system('python ~/github/HalphaImaging/python3/uat_astr_mosaic.py --s --filestring ghz --psfex')
        # run psfex on convolved images

        os.system('ls ghz*o00.cat > gauss_list')
        os.system('psfex @gauss_list -XML_NAME gpsfex.xml')
        #

        gimage_names, gimage_fwhm = get_fwhm_psfex(psfout='gpsfex.xml')


        print('#######################################')
        print('COMPARING ORIGINAL AND FINAL FWHM')
        print('#######################################')
        print("image   fwhm    gfwhm")
        print("------------------------")
        for i in range(len(image_names)):
            print(f"{image_names[i]}: {image_fwhm[i]:.2f} {gimage_fwhm[i]:.2f} ratio={gimage_fwhm[i]/fwhm_max:.3f}")

            # write to a file
        newtab = Table([image_names,image_fwhm,gimage_fwhm])
        newtab.write("psfex_original_and_convolved_fwhm.csv",format='csv',overwrite=True)
