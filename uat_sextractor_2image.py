#!/usr/bin/env python

'''
GOAL:
  The goal of this program is to run SExtractor in two image mode, to detect objects on the R-band image and
  apply the same apertures to the Halpha image.

PROCEDURE:
  - copy default setup files to current directory
  - Run sextractor on each image

EXAMPLE:
   In the directory containing all flattened objects with fixed headers to run sextractor type in the command line:
      '/Users/alfalfa/Github/HalphaImaging/uat_sextractor_2image.py --s'(or whatever the path is to where this program is stored)

   To get swarp to create aligned images in multiple bans (e.g. Halpha and R-band), do the following
    uat_astr_mosaic.py --swarp --l A1367-h02_ha12

    uat_astr_mosaic.py --swarp --l A1367-h02_R --refimage 'A1367-h02_ha12.coadd.fits'

    uat_astr_mosaic.py --swarp --l A1367-h02_ha12 --refimage 'A1367-h02_ha12.coadd.fits'


WHAT THIS CODE DOES:
INPUT/OUPUT:
REQUIRED MODULES:
EXTRA NOTES:
WRITTEN BY:
Rose Finn, 04 Jan 2017

'''
import glob
import os
from astropy.io import fits
from astropy.wcs import WCS
import argparse
import subprocess
from matplotlib import pyplot as plt
import numpy as np

parser = argparse.ArgumentParser(description ='Run sextractor in two-image mode')
parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Run sextractor to create object catalogs')
parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of default config files')
parser.add_argument('--image1',dest = 'image1', default = None,  help = 'image used to define apertures')
parser.add_argument('--image2',dest = 'image2', default = None,  help = 'image used to for measuring phot based on image1 (typically this is the Halpha image)')
parser.add_argument('--plot',dest = 'plot', default = False, action = 'store_true', help = 'make diagnostic plots')
parser.add_argument('--imagedir',dest = 'imagedir', default = '/Users/rfinn/research/HalphaGroups/plots/', help = 'directory for saving plots')


args = parser.parse_args()

# get input files
#print 'cp ' +args.d + '/default.* .'
os.system('cp ' +args.d + '/default.* .')
#files = sorted(glob.glob(args.filestring))

#nfiles = len(files)
i = 1


if args.s:

    print 'RUNNING SEXTRACTOR'
    t = args.image1.split('.fits')
    froot1 = t[0]
    os.system('sex ' + args.image1+','+args.image1 + ' -c default.sex.hdi -CATALOG_NAME ' + froot1 + '.cat')
    os.rename('check.fits', froot1 + 'check.fits')
    # run on second image
    t = args.image2.split('.fits')
    froot2 = t[0]
    os.system('sex ' + args.image1+','+args.image2 + ' -c default.sex.hdi -CATALOG_NAME ' + froot2 + '.cat')
    os.rename('check.fits', froot2 + 'check.fits')

    #plt.figure()
    #plt.imshow(froot1+'check.fits')
    #plt.title('Image 1')
    #plt.figure()
    #plt.imshow(froot2+'check.fits')
    #plt.title('Image 2')

if args.plot:
    t = args.image1.split('.fits')
    froot1 = t[0]
    cat1 = fits.getdata(froot1+'.cat',2)
    t = args.image2.split('.fits')
    froot2 = t[0]
    cat2 = fits.getdata(froot2+'.cat',2)
    plt.figure(figsize=(6,4))
    plt.subplots_adjust(bottom=.2,left=.15,right=.95)
    x = cat2.FLUX_AUTO
    y = cat2.FLUX_AUTO/cat1.FLUX_AUTO
    flag = cat1.FLAGS == 0
    plt.plot(x[flag],y[flag],'k.')
    #plt.axis([0,500,.02,.08])
    x1,x2 = plt.xlim()
    plt.xlim(.2,x2)
    ave = np.median(y[(x> 50) & flag])
    std = np.std(y[(x > 50) & flag])
    plt.axhline(y=ave)
    print '%.4f (%.4f)'%(ave,std)
    plt.ylabel('$Flux (Halpha)/Flux(R) $',fontsize=20)
    plt.xlabel('$Flux(R) \ (ADU)$',fontsize=20)
    plt.text(20,.07,'$ ratio = %.4f (%.4f)$'%(ave,std),fontsize=16)
    plt.gca().set_xscale('log')
    t = args.image2.split('.coadd')
    plt.title(t[0],fontsize=20)
    plt.savefig(args.imagedir+t[0]+'-filter-ratio.png')
