#!/usr/bin/env python


'''

GOAL:

PROCEDURE:
   scale R band by a factor
   subtract from Halpha
   save resultant imaage
  
EXAMPLE:

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

args = parser.parse_args()


r = fits.getdata(args.r)
ha = fits.getdata(args.ha)
scale = float(args.scale)

adjust_scale = True
while adjust_scale:
    cs = ha  - scale*r
    v1,v2=scoreatpercentile(ha,[5.,95.])
    plt.figure(1,figsize=(10,4))
    plt.clf()
    plt.subplots_adjust(hspace=0,wspace=0)
    plt.subplot(1,3,1)
    plt.imshow(r,cmap='gray_r',vmin=v1,vmax=v2)
    plt.title('R')
    plt.subplot(1,3,2)
    plt.imshow(ha,cmap='gray_r',vmin=v1,vmax=v2)
    plt.title('Halpha')
    plt.gca().set_yticks(())
    plt.subplot(1,3,3)
    plt.imshow(cs,cmap='gray_r',vmin=v1,vmax=v2)
    plt.gca().set_yticks(())
    plt.title('contsub, scale = %4.3f'%(scale))
    plt.draw()
    plt.show()
    
    
    t=raw_input('enter new value for scale factor;\n w to write output and quit;\n q to quit without saving\n')
    try:
        scale=float(t)
    except ValueError:
        adjust_scale = False
        if t.find('q') > -1:
            break
        if t.find('w') > -1:
            outimage = args.ha.split('-Ha')[0]+'-cs.fits'
            newfile = fits.PrimaryHDU()
            newfile.data = cs
            newfile.header = ha.header
            fits.writeto(outimage, newfile.data, header = newfile.header, clobber=True)
