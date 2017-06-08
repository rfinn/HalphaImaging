#!/usr/bin/env python


'''

GOAL:

To create a continuum-subtracted image given
(1) an R-band and Halpha image, and
(2) the scale factor to apply to the R-band image.

PROCEDURE:

- scale R band by a factor
- subtract from Halpha
- iterate until the you are happy
- save resultant continuum-subtracted image if desired
  
EXAMPLES:

(1) to run within ipython type:

%run ~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale .055 --mosaic

The mosaic flag increases the size of the images that are displayed so you can see the results better.

(2) to run from the command line type:

~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale .055

(3) to run on cutout images rather than mosaics:

 ~/github/HalphaImaging/uat_subtract_continuum.py --cluster A1367 --scale 0.044 --id 113364

INPUT/OUPUT:

REQUIRED MODULES:

astropy
argparse
matplotlib
scipy

EXTRA NOTES:

WRITTEN BY:
   Rose Finn

'''
from astropy.io import fits
import argparse
from argparse import RawDescriptionHelpFormatter
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

parser = argparse.ArgumentParser(description ='This program subtracts scaled R-band image from Halpha.\n \nTo subtract mosaics:\n~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale 0.0445 --mosaic \n\nTo subtract cutouts:\n~/github/HalphaImaging/uat_subtract_continuum.py --cluster A1367 --scale 0.044 --id 113364', formatter_class=RawDescriptionHelpFormatter)
parser.add_argument('--cluster',dest = 'cluster', default = None, help = 'Cluster prefix of image for continuum subtraction.  Use this if you are subtracting continuum from a cutout rather than a mosaic.')
parser.add_argument('--id',dest = 'id', default = None, help = 'NSAID of image for continuum subtraction.  Use this if you are subtracting continuum from a cutout rather than a mosaic.')
parser.add_argument('--r', dest = 'r', default = None, help = 'R-band image.  Use this if you are subtracting mosaic images rather than cutouts.')
parser.add_argument('--ha', dest = 'ha', default = None, help = 'Halpha image.  Use this if you are subtracting mosaic images rather than cutouts.')
parser.add_argument('--scale', dest = 'scale', default = 0.0445, help = 'factor to scale R-band image by before subtracting from Halpha image')
parser.add_argument('--mosaic', dest = 'mos', default = False, action = 'store_true', help = 'set this if subtracting mosaic images rather than cutouts.  It will make the figure bigger so you can see the mosaics better.')
args = parser.parse_args()

figure_size=(10,4)
if args.mos:
    figure_size=(14,6)
    
if args.id:
    rimage = args.cluster+'-'+args.id+'-R.fits'
    haimage = args.cluster+'-'+args.id+'-Ha.fits'
    id = args.id
else:
    rimage = args.r
    haimage = args.ha
    t = args.r.split('-')
    id = t[0]
r,r_header = fits.getdata(rimage,header=True)
ha,ha_header = fits.getdata(haimage,header=True)

scale = float(args.scale)

adjust_scale = True
while adjust_scale:
    plt.close('all')
    cs = ha  - scale*r
    v1,v2=scoreatpercentile(ha,[.5,99.])
    
    plt.figure(1,figsize=figure_size)
    plt.clf()
    plt.subplots_adjust(hspace=0,wspace=0)
    #Halpha plus continuum
    plt.subplot(1,3,1)
    plt.imshow(ha,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
    plt.title('Halpha + cont')
    #plt.gca().set_yticks(())
    plt.xlabel('NSA ID '+id,fontsize=14)
    #R
    plt.subplot(1,3,2)
    plt.imshow(r,cmap='gray_r',vmin=v1/scale,vmax=v2/scale,origin='lower')
    plt.title('R')
    plt.gca().set_yticks(())
    #Continuum subtracted image
    plt.subplot(1,3,3)
    v3,v4=scoreatpercentile(cs,[1.,99.])
    plt.imshow(cs,origin='lower',cmap='gray_r',vmin=v1,vmax=v2)
    plt.gca().set_yticks(())
    plt.title('contsub, scale = %4.4f'%(scale))
    plt.show(block=False)
    
    
    t=raw_input('enter new value for scale factor;\n   w to write output and quit;\n   q to quit without saving\n')
    try:
        scale=float(t)
    except ValueError:
        adjust_scale = False
        if t.find('q') > -1:
            break
        if t.find('w') > -1:
            if haimage.find('Ha') > -1:
                outimage = haimage.split('-Ha')[0]+'-CS.fits'
            elif haimage.find('ha') > -1:
                outimage = haimage.split('_ha')[0]+'-CS.fits'
            newfile = fits.PrimaryHDU()
            newfile.data = cs
            newfile.header = ha_header
            fits.writeto(outimage, newfile.data, header = newfile.header, clobber=True)
            output = outimage.split('.fits')
            plt.savefig(output[0]+'.png')
