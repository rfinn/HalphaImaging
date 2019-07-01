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

~/github/HalphaImaging/uat_make_cutouts.py --image A1367-h02_ha12.coadd.fits --catalog ~/research/NSA/nsa_v0_1_2.fits --filter Ha --nhalpha 12 

~/github/HalphaImaging/uat_make_cutouts.py --image A1367-h02_R.coadd.fits --catalog ~/research/NSA/nsa_v0_1_2.fits --filter R --nhalpha 12 

~/github/HalphaImaging/uat_make_cutouts.py --image A1367-h02-CS.fits --catalog ~/research/NSA/nsa_v0_1_2.fits --filter CS --nhalpha 12 


To run on Virgo filament imaging:

~/github/HalphaImaging/uat_make_cutouts.py --image pointing-13_R.coadd.fits --catalog ~/github/Virgo/tables/nsa.virgo.fits --filter R --nhalpha 4 --prefix pointing13

~/github/HalphaImaging/uat_make_cutouts.py --image pointing-13_ha4.coadd.fits --catalog ~/github/Virgo/tables/nsa.virgo.fits --filter ha4 --nhalpha 4 --prefix pointing13 --plot

    

EXAMPLE:

    uat_make_cutouts.py --image A1367-h02_R.coadd.fits --filter R --nhalpha 12 --cluster coma

    uat_make_cutouts.py --image A1367-h02_R.coadd.fits --catalog /Users/rfinn/research/NSA/nsa_v0_1_2.fits --filter R --nhalpha 12 --plot

NATASHA NOTES:
dest is the user parameters and the argument is just used in the code

REQUIRED MODULES:

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



parser = argparse.ArgumentParser(description ='Get cutouts for NSA galaxies within field of view of mosaic and redshift range of designated H-alpha filter')
parser.add_argument('--image', dest = 'image', default = None, help = 'mosaic/HDI image to make cutouts from')
parser.add_argument('--catalog', dest = 'catalog', default = '/home/share/catalogs/nsa_v0_1_2.fits', help = 'full path to the NSA catalog')
parser.add_argument('--Haimage', dest = 'Haimage', default = None, help = 'H-alpha HDI/mosaic image to make cutouts from')
parser.add_argument('--filter',dest = 'filter', default ='R', help = 'Filter for the input mosaic image (e.g. r, R, Ha).  Default value is R.')
parser.add_argument('--nhalpha',dest = 'nhalpha', default ='12', help = 'H-alpha filter number (e.g. 4, 8, 12 or 16).  Default value is 12.')
parser.add_argument('--Rscale',dest = 'scale', default =15., help = 'cutout size = (scale x Re, scale x Re) - increase scale to increase size of cutout.  Default value is 15.')
parser.add_argument('--prefix',dest = 'prefix', default ='A1367', help = 'cluster name to preprend to cutout image names (no spaces!).  Default value is A1367.')
parser.add_argument('--plot', dest = 'plot', default = False, action = 'store_true', help = 'plot cutouts and position wrt mosaic.  Default value is False.')
#parser.add_argument('--both',dest = 'both', default = False, help = 'Runs cutouts for both Ha and R images at same time')
#parse.add_argument('--filter2',dest = 'filter2',default = 'Ha', help = 'Filter for the mosaic image of Ha, helps to run both at same time')
#parser.add_argument('--l', dest = 'l', default = False, help = 'List of images to input to swarp')

args = parser.parse_args()



# setting up filter information
# dictionary of Halpha filters
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.,'INT197':6520.5}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.,'INT197':6615.5}

Zmax=(((lmax[args.nhalpha])/6563.)-1)
Zmin=(((lmin[args.nhalpha])/6563.)-1)
print 'Galaxies detectable in Halpha have redshifts between ',Zmin,' and ', Zmax

def makecuts(image,imagefilter):
    catdat= fits.getdata(args.catalog)
    print 'Cutting out', image    
    
    zFlag = (catdat.Z > Zmin) & (catdat.Z < Zmax)
    
    f = fits.open(image)
    prihdr = f[0].header
    n2,n1 = f[0].data.shape
    
    w= WCS(image)
    px,py = w.wcs_world2pix(catdat.RA,catdat.DEC,1)
    onimageflag=(px < n1) & (px >0) & (py < n2) & (py > 0)
    
    keepflag=zFlag & onimageflag
    RA=catdat.RA[keepflag]
    DEC=catdat.DEC[keepflag]
    radius=catdat.SERSIC_TH50[keepflag]
    IDNUMBER=catdat.NSAID[keepflag]
    redshift = catdat.Z[keepflag]
    zdist = catdat.ZDIST[keepflag]

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
            cutout = Cutout2D(f[0].data, position, size, wcs=w, mode='trim') #require entire image to be on parent image
        except astropy.nddata.utils.PartialOverlapError:# PartialOverlapError:
            print 'galaxy is only partially covered by mosaic - skipping ',IDNUMBER[i]
            continue
        except astropy.nddata.utils.NoOverlapError:# PartialOverlapError:
            print 'galaxy is not covered by mosaic - skipping ',IDNUMBER[i]
            continue
        if args.plot:
            plt.figure()
            plt.imshow(f[0].data, origin='lower',cmap='gray', norm=LogNorm())
            cutout.plot_on_original(color='white')
            plt.show()
            r = raw_input('type any key to continue (p to skip plotting) \n')
            if r.find('p') > -1:
                args.plot = False
        # figure out how to save the cutout as fits image
        ((ymin,ymax),(xmin,xmax)) = cutout.bbox_original
        outimage = args.prefix+'-'+(str(IDNUMBER[i])+'-'+ args.filter+".fits")
        newfile = fits.PrimaryHDU()
        newfile.data = f[0].data[ymin:ymax,xmin:xmax]
        newfile.header = f[0].header
        newfile.header.update(w[ymin:ymax,xmin:xmax].to_header())
        newfile.header.set('REDSHIFT',float('{:.6f}'.format(redshift[i])))
        newfile.header.set('ZDIST',float('{:.6f}'.format(zdist[i])))
        newfile.header.set('NSAID',float('{:d}'.format(IDNUMBER[i])))
        
        fits.writeto(outimage, newfile.data, header = newfile.header, overwrite=True)
    return cutout
if __name__ == '__main__':
    imcutout = makecuts(args.image,args.filter)

