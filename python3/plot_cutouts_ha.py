#!/usr/bin/env python

from astropy.io import fits
import argparse
from argparse import RawDescriptionHelpFormatter
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
import sys

vmin = .5
vmax = 99.5

class cutouts():
    def __init__(self, rimage):
        self.r_name = rimage
        if self.r_name.find('_R') > -1:
            split_string = '_R.fits'

        if self.r_name.find('-R') > -1:
            split_string = '-R.fits'
        elif self.r_name.find('-r') > -1:
            split_string = '_r.fits'

        self.rootname = self.r_name.split(split_string)[0]            
        self.r = fits.getdata(self.r_name)
        self.ha = fits.getdata(self.rootname+'-Ha.fits')
        self.cs = fits.getdata(self.rootname+'-CS.fits')
    def plotcutouts(self):
        figure_size=(10,4)
        
        v1,v2=scoreatpercentile(self.ha,[vmin,vmax])#.5,99
    
        plt.figure(figsize=figure_size)
        plt.clf()
        plt.subplots_adjust(hspace=0,wspace=0)
        #Halpha plus continuum
        plt.subplot(1,3,1)
        plt.imshow(self.ha,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
        plt.title(r'$H\alpha + cont$',fontsize=14)
        #plt.gca().set_yticks(())

        #R
        plt.subplot(1,3,2)
        v1,v2=scoreatpercentile(self.r,[vmin,vmax])#.5,99        
        plt.imshow(self.r,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
        plt.title(r'$R$',fontsize=14)
        plt.gca().set_yticks(())
        #Continuum subtracted image
        plt.xlabel(self.rootname,fontsize=14)        
        plt.subplot(1,3,3)
        v1,v2=scoreatpercentile(self.cs,[vmin,vmax])
        plt.imshow(self.cs,origin='lower',cmap='gray_r',vmin=v1,vmax=v2)
        plt.gca().set_yticks(())
        plt.title(r'$H\alpha$',fontsize=14)
        plt.show(block=False)
        plt.savefig(self.rootname+'-cutouts.png')
        plt.savefig(self.rootname+'-cutouts.pdf')        


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='This program subtracts scaled R-band image from Halpha.\n \nTo subtract mosaics:\n~/github/HalphaImaging/uat_subtract_continuum.py --r A1367-h02_R.coadd.fits --ha A1367-h02_ha12.coadd.fits --scale 0.0445 --mosaic \n\nTo subtract cutouts:\n~/github/HalphaImaging/uat_subtract_continuum.py --cluster A1367 --scale 0.044 --id 113364', formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--r', dest = 'r', default = None, help = 'R-band image.  This is all you need to provide if images are named: blah_R.fits, blah_CS.fits, blah_Ha.fits')
    parser.add_argument('--ha', dest = 'ha', default = None, help = 'Halpha image.  Use this if you are subtracting mosaic images rather than cutouts.')
    args = parser.parse_args()

    if args.r is None:
        print('you must supply the r-band image name!')
        print('try again')
        sys.exit()
    c = cutouts(args.r)
    c.plotcutouts()
