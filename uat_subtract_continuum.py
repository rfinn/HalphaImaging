#!/usr/bin/env python


'''

GOAL:

PROCEDURE:
   scale R band by a factor
   subtract from Halpha
   save resultant imaage
  
EXAMPLE:
   from within ipython type:

   %run ~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale .055

INPUT/OUPUT:

REQUIRED MODULES:

EXTRA NOTES:

WRITTEN BY:
   Rose Finn

'''
from astropy.io import fits
import argparse
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

parser = argparse.ArgumentParser(description ='subtract scaled R-band image from Halpha')
parser.add_argument('--r', dest = 'r', default = None, help = 'R-band image')
parser.add_argument('--ha', dest = 'ha', default = None, help = 'R-band image')
parser.add_argument('--scale', dest = 'scale', default = 0.06, help = 'R-band image')
parser.add_argument('--mosaic', dest = 'mos', default = False, action = 'store_true', help = 'specifies mosaic image')
args = parser.parse_args()

figure_size=(10,4)
if args.mos:
    figure_size=(15,8)
    

r,r_header = fits.getdata(args.r,header=True)
ha,ha_header = fits.getdata(args.ha,header=True)
scale = float(args.scale)

adjust_scale = True
while adjust_scale:
    cs = ha  - scale*r
    v1,v2=scoreatpercentile(ha,[5.,95.])
    plt.figure(1,figsize=figure_size)
    plt.clf()
    plt.subplots_adjust(hspace=0,wspace=0)
    #Halpha plus continuum
    plt.subplot(1,3,1)
    plt.imshow(ha,cmap='gray_r',vmin=v1,vmax=v2)
    plt.title('Halpha + cont')
    plt.gca().set_yticks(())
    plt.gca().invert_yaxis(())
    #R
    plt.subplot(1,3,2)
    plt.imshow(r,cmap='gray_r',vmin=v1,vmax=v2)
    plt.title('R')
    plt.gca().set_yticks(())
    plt.gca().invert_yaxis((()
    #Continuum subtracted image
    plt.subplot(1,3,3)
    plt.imshow(cs,cmap='gray_r',vmin=v1,vmax=v2)
    plt.title('contsub, scale = %4.3f'%(scale))
    plt.gca().set_yticks(())
    plt.gca().invert_yaxis(())
    plt.draw()
    plt.show()
    
    
    t=raw_input('enter new value for scale factor;\n   w to write output and quit;\n   q to quit without saving\n')
    try:
        scale=float(t)
    except ValueError:
        adjust_scale = False
        if t.find('q') > -1:
            break
        if t.find('w') > -1:
            outimage = args.ha.split('-Ha')[0]+'-CS.fits'
            newfile = fits.PrimaryHDU()
            newfile.data = cs
            newfile.header = ha_header
            fits.writeto(outimage, newfile.data, header = newfile.header, clobber=True)
