#!/usr/bin/env python

from astropy.io import fits
import argparse
from argparse import RawDescriptionHelpFormatter
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile


vmin = .5
vmax = 98

class cutouts():
    def __init__(self, rimage):
        self.r = rimage
        if self.r.find('_R') > -1:
            self.rootname = self.r.split('_R.fits')
        elif self.r.find('-r') > -1:
            self.rootname = self.r.split('_r.fits')            
        self.ha = self.rootname+'_Ha.fits'
        self.cs = self.rootname+'_CS.fits'        
    def plotcutouts(self):
        figure_size=(10,4)
        
        v1,v2=scoreatpercentile(self.ha,[vmin,vmax])#.5,99
    
        plt.figure(1,figsize=figure_size)
        plt.clf()
        plt.subplots_adjust(hspace=0,wspace=0)
        #Halpha plus continuum
        plt.subplot(1,3,1)
        plt.imshow(ha,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
        plt.title('Halpha + cont')
        #plt.gca().set_yticks(())
        plt.xlabel(self.prefix,fontsize=14)
        #R
        plt.subplot(1,3,2)
        v1,v2=scoreatpercentile(self.r,[vmin,vmax])#.5,99        
        plt.imshow(r,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
        plt.title('R')
        plt.gca().set_yticks(())
        #Continuum subtracted image
        plt.subplot(1,3,3)
        v1,v2=scoreatpercentile(self.cs,[vmin,vmax])
        plt.imshow(cs,origin='lower',cmap='gray_r',vmin=v1,vmax=v2)
        plt.gca().set_yticks(())
        plt.title('contsub')
        plt.show(block=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='This program subtracts scaled R-band image from Halpha.\n \nTo subtract mosaics:\n~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale 0.0445 --mosaic \n\nTo subtract cutouts:\n~/github/HalphaImaging/uat_subtract_continuum.py --cluster A1367 --scale 0.044 --id 113364', formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--r', dest = 'r', default = None, help = 'R-band image.  This is all you need to provide if images are named: blah_R.fits, blah_CS.fits, blah_Ha.fits')
    parser.add_argument('--ha', dest = 'ha', default = None, help = 'Halpha image.  Use this if you are subtracting mosaic images rather than cutouts.')
    args = parser.parse_args()

    if args.r is None:
        print('you must supply the r-band image name!')
        print('try again')
    else:
        cimage = cutouts(args.r)
