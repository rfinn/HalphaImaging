#! /usr/bin/env python 

"""

GOAL:

    to create postage-stamp cutouts of galaxies

PROCEDURE:

    read catalog containing RA, Dec and some measure of galaxy size (e.g. NSA THETA_50)
    keep RA, DEC for galaxies that are on image and in the right redshift slice so that Halpha is viewable through specified filter
    feed results RA, DEC and Size to Cutout2D

INPUTS:

    catalog - catalog or subcatalog of NASA Sloan atlas
    halpha filter number - this is the filter number used to refer to the Halpha filter for your group/cluster (e.g. Ha8 would be 8)

USAGE: 

When running from macbook (not coma):

~/github/HalphaImaging/uat_make_cutouts.py --image A1367-h02_ha12.coadd.fits --catalog ~/research/NSA/nsa_v0_1_2.fits --filter Ha --nhalpha 12 --cluster A1367

To run on Virgo filament imaging:

~/github/HalphaImaging/uat_make_cutouts.py --image pointing-13_R.coadd.fits --catalog ~/github/Virgo/tables/nsa.virgo.fits --filter R --nhalpha 4 --prefix pointing13

~/github/HalphaImaging/uat_make_cutouts.py --image pointing-13_ha4.coadd.fits --catalog ~/github/Virgo/tables/nsa.virgo.fits --filter ha4 --nhalpha 4 --prefix pointing13 --plot

%run ~/github/HalphaImaging/uat_make_both_cutouts.py --Rimage pointing-4_r.coadd.fits --Haimage pointing-4_ha4.coadd.fits --catalog ~/github/Virgo/tables/nsa.virgo.fits --filter1 R --filter2 Ha --nhalpha 4 --prefix pointing4 

REQUIRED MODULES:

astropy
matplotlib
numpy

NOTES:
    You may need to change the variables Haimage, Rimage and Hawcimage if your image names are different

REFERENCES:

http://docs.astropy.org/en/stable/nddata/utils.html

INT filters from here
http://catserver.ing.iac.es/filter/list.php?instrument=WFC

Details for INT filter # 197
central wavelength = 6568
width = 95
lmin = 6568 - .5*95 = 6540.5
lmax = 6568 + 0.5*95 = 6615.5



"""
import numpy as np
from matplotlib import pyplot as plt
import astropy
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, FK5
import astropy.units as u
from astropy.nddata.utils import Cutout2D
import argparse
from matplotlib.colors import LogNorm
#try:
#    import montage
#    import os
#    CanUseMontage=True
#except ImportError:
#    CanUseMontage=False
#except Exception:
#    CanUseMontage=False



parser = argparse.ArgumentParser(description ='Get cutouts for NSA galaxies within field of view of mosaic and redshift range of designated H-alpha filter \n example: \n %run ~/github/HalphaImaging/uat_make_both_cutouts.py --Rimage pointing-4_r.coadd.fits --Haimage pointing-4_ha4.coadd.fits --catalog ~/github/Virgo/tables/nsa.virgo.fits --filter1 R --filter2 Ha --nhalpha 4 --prefix pointing4  ')
parser.add_argument('--Rimage', dest = 'Rimage', default = None, help = 'R-band HDI/mosaic image to make cutouts from')
parser.add_argument('--Haimage', dest = 'Haimage', default = None, help = 'H-alpha HDI/mosaic image to make cutouts from')
parser.add_argument('--catalog', dest = 'catalog', default = '/home/share/catalogs/nsa_v0_1_2.fits', help = 'full path to the NSA catalog')
parser.add_argument('--filter1',dest = 'filter1', default ='R', help = 'Filter for the input mosaic image (e.g. r, R, Ha).  Default value is R.')
parser.add_argument('--filter2',dest = 'filter2', default ='Ha', help = 'Filter for the input mosaic image (e.g. r, R, Ha).  Default value is Ha.')
parser.add_argument('--nhalpha',dest = 'nhalpha', default ='4', help = 'H-alpha filter number (e.g. 4, 8, 12, 16 or INT197).  Default value is 4.')
parser.add_argument('--Rscale',dest = 'scale', default =16., help = 'cutout size = (scale x Re, scale x Re) - increase scale to increase size of cutout.  Default value is 16.')
parser.add_argument('--prefix',dest = 'prefix', default ='A1367', help = 'cluster name to preprend to cutout image names (no spaces!).  Default value is A1367.')
parser.add_argument('--plot', dest = 'plot', default = False, action = 'store_true', help = "add '--plot ' to plot cutouts and position wrt mosaic.  Default value is False.")


args = parser.parse_args()


# setting up filter information
#dictionary of Halpha filters
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.,'INT197':6540.5}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.,'INT197':6615.5}

# convert filter min and max wavelength to redshift for halpha emission
Zmax=(((lmax[args.nhalpha])/6563.)-1)
Zmin=(((lmin[args.nhalpha])/6563.)-1)
print 'Galaxies detectable in Halpha have redshifts between ',Zmin,' and ', Zmax


def makebothcuts(Rimage,filter1,Haimage,filter2):
    catdat= fits.getdata(args.catalog)
    print 'Cutting out', Rimage
    print 'Cutting out', Haimage
   
    
    zFlag = (catdat.Z > Zmin) & (catdat.Z < Zmax)
    
    f = fits.open(Rimage)
    g = fits.open(Haimage)
    prihdr = f[0].header
    prihdr1 = g[0].header
    n2,n1 = f[0].data.shape #should be same for Ha too, maybe? IDK
    n4,n3 = g[0].data.shape 
    
    w= WCS(Rimage)#OF R IMAGE, SO THAT HA MATCHES WCS OF R, SO THEY'RE THE SAME
    px,py = w.wcs_world2pix(catdat.RA,catdat.DEC,1)
    onimageflag=(px < n1) & (px >0) & (py < n2) & (py > 0)
    
    keepflag=zFlag & onimageflag
    RA=catdat.RA[keepflag]
    DEC=catdat.DEC[keepflag]
    radius=catdat.SERSIC_TH50[keepflag]
    IDNUMBER=catdat.NSAID[keepflag]
    print 'number of galaxies to keep = ', sum(keepflag)

#    if args.region_file:
        
    for i in range(len(RA)):

        if (radius[i]<.01):
            size=120.
        else:
            size=float(args.scale)*radius[i]
            
        position = SkyCoord(ra=RA[i],dec=DEC[i],unit='deg')
        size = u.Quantity((size, size), u.arcsec)
        #print image, radius[i], position, size
        #cutout = Cutout2D(fdulist[0].data, position, size, wcs=w, mode='strict') #require entire image to be on parent image
        try:
            cutoutR = Cutout2D(f[0].data, position, size, wcs=w, mode='trim') #require entire image to be on parent image
            cutoutHa = Cutout2D(g[0].data, position, size, wcs=w, mode = 'trim')
        except astropy.nddata.utils.PartialOverlapError:# PartialOverlapError:
            print 'galaxy is only partially covered by mosaic - skipping ',IDNUMBER[i]
            continue
        except astropy.nddata.utils.NoOverlapError:# PartialOverlapError:
            print 'galaxy is not covered by mosaic - skipping ',IDNUMBER[i]
            continue
        if args.plot:
            plt.figure() #R cutout
            plt.imshow(f[0].data, origin='lower',cmap='gray', norm=LogNorm())
            cutoutR.plot_on_original(color='white')
            plt.show()
            r = raw_input('type any key to continue (p to skip plotting) \n')
            if r.find('p') > -1:
                args.plot = False
            plt.figure() #Ha cutout
            plt.imshow(g[0].data, origin='lower',cmap='gray', norm=LogNorm())
            cutoutHa.plot_on_original(color='white')
            plt.show()
            r = raw_input('type any key to continue (p to skip plotting) \n')
            if r.find('p') > -1:
                args.plot = False
        # saving R Cutout as fits image
        ((ymin,ymax),(xmin,xmax)) = cutoutR.bbox_original
        outimage = args.prefix+'-'+(str(IDNUMBER[i])+'-'+ args.filter1+".fits")
        newfile = fits.PrimaryHDU()
        newfile.data = f[0].data[ymin:ymax,xmin:xmax]
        newfile.header = f[0].header
        newfile.header.update(w[ymin:ymax,xmin:xmax].to_header())
        
        fits.writeto(outimage, newfile.data, header = newfile.header, clobber=True)
        # saving Ha Cutout as fits image
        ((ymin1,ymax1),(xmin1,xmax1)) = cutoutHa.bbox_original
        outimage1 = args.prefix+'-'+(str(IDNUMBER[i])+'-'+ args.filter2+".fits")
        newfile1 = fits.PrimaryHDU()
        newfile1.data = g[0].data[ymin1:ymax1,xmin1:xmax1]
        newfile1.header = g[0].header
        newfile1.header.update(w[ymin1:ymax1,xmin1:xmax1].to_header())
        newfile.header.set('REDSHIFT',float('{:.6f}'.format(redshift[i])))
        newfile.header.set('ZDIST',float('{:.6f}'.format(zdist[i])))
        newfile.header.set('NSAID',float('{:d}'.format(IDNUMBER[i])))
        

        fits.writeto(outimage1, newfile1.data, header = newfile1.header, clobber=True)
        if args.plot:
           plt.figure()
           plt.imshow(g[0].data, origin='lower',cmap='gray', norm=LogNorm())
           cutoutHa.plot_on_original(color='white')
           plt.show()
           r = raw_input('type any key to continue (p to skip plotting) \n')
           if r.find('p') > -1:
                args.plot = False
        # figure out how to save the cutout as fits image
        #((ymin,ymax),(xmin,xmax)) = cutoutHa.bbox_original
       # outimage = args.prefix+'-'+(str(IDNUMBER[i])+'-'+ args.filter2+".fits")
       # newfile = fits.PrimaryHDU()
       # newfile.data = g[0].data[ymin:ymax,xmin:xmax]
       # newfile.header = g[0].header
       # newfile.header.update(w[ymin:ymax,xmin:xmax].to_header())
        
       # fits.writeto(outimage, newfile.data, header = newfile.header, clobber=True)
if __name__ == '__main__':
    imcutout1 = makebothcuts(args.Rimage,args.filter1,args.Haimage,args.filter2)
    #imcutout2 = makebothcuts(args.Rimage,args.filter1,args.Haimage,args.filter2)

