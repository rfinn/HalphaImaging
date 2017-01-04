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

~/github/HalphaImaging/uat_make_cutouts.py --image A1367-h02_R.coadd.fits --catalog ~/research/NSA/nsa_v0_1_2.fits --filter R --nhalpha 12 --cluster A1367

~/github/HalphaImaging/uat_make_cutouts.py --image A1367-h02-CS.fits --catalog ~/research/NSA/nsa_v0_1_2.fits --filter CS --nhalpha 12 --cluster A1367

If running on coma, you don't need to specify the catalog path - the default should be the correct value.

    uat_make_cutouts.py catalog image Halpha_filter_number


EXAMPLE:

    uat_make_cutouts.py --image A1367-h02_R.coadd.fits --filter R --nhalpha 12 --cluster coma

    uat_make_cutouts.py --image A1367-h02_R.coadd.fits --catalog /Users/rfinn/research/NSA/nsa_v0_1_2.fits --filter R --nhalpha 12 --plot


REQUIRED MODULES:

NOTES:
    You may need to change the variables Haimage, Rimage and Hawcimage if your image names are different

REFERENCES:

http://docs.astropy.org/en/stable/nddata/utils.html

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
parser.add_argument('--filter',dest = 'filter', default ='R', help = 'Filter for the input mosaic image (e.g. r, R, Ha).  Default value is R.')
parser.add_argument('--nhalpha',dest = 'nhalpha', default ='12', help = 'H-alpha filter number (e.g. 4, 8, 12 or 16).  Default value is 12.')
parser.add_argument('--scale',dest = 'scale', default =15., help = 'cutout size = (scale x Re, scale x Re) - increase scale to increase size of cutout.  Default value is 15.')
parser.add_argument('--cluster',dest = 'cluster', default ='A1367', help = 'cluster name to preprend to cutout image names (no spaces!).  Default value is A1367.')
parser.add_argument('--plot', dest = 'plot', default = False, action = 'store_true', help = 'plot cutouts and position wrt mosaic.  Default value is False.')
#parser.add_argument('--l', dest = 'l', default = False, help = 'List of images to input to swarp')

args = parser.parse_args()



# setting up filter information
#dictionary of Halpha filters
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.}

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
        outimage = args.cluster+'-'+(str(IDNUMBER[i])+'-'+ args.filter+".fits")
        newfile = fits.PrimaryHDU()
        newfile.data = f[0].data[ymin:ymax,xmin:xmax]
        newfile.header = f[0].header
        newfile.header.update(w[ymin:ymax,xmin:xmax].to_header())
        
        fits.writeto(outimage, newfile.data, header = newfile.header, clobber=True)
    return cutout
if __name__ == '__main__':
    imcutout = makecuts(args.image,args.filter)


########## OLD STUFF ############

#Name of images
Haimage=('Hacs_final.fits')
Rimage= ('R_final.fits')
Hawcimage=('Ha_final.fits')
#Hawc is image with continuum 


def get_cutout(image,RA,Dec,size):
    print 'getting cutouts!'


def cutout(filename, xc, yc, xw=25, yw=25, units='pixels', outfile=None,
        clobber=True, useMontage=False, coordsys='celestial', verbose=False):
    
    #Inputs:
        #file  - .fits filename or pyfits HDUList (must be 2D)
        #xc,yc - x and y coordinates in the fits files' coordinate system (CTYPE)
        #xw,yw - x and y width (pixels or wcs)
        #units - specify units to use: either pixels or wcs
        #outfile - optional output file
    

    if isinstance(filename,str):
        file = pyfits.open(filename)
        opened=True
    elif isinstance(filename,pyfits.HDUList):
        file = filename
        opened=False
    else:
        raise Exception("cutout: Input file is wrong type (string or HDUList are acceptable).")

    head = file[0].header.copy()

    if head['NAXIS'] > 2:
        raise DimensionError("Too many (%i) dimensions!" % head['NAXIS'])
    cd1 = head.get('CDELT1') if head.get('CDELT1') else head.get('CD1_1')
    cd2 = head.get('CDELT2') if head.get('CDELT2') else head.get('CD2_2')
    if cd1 is None or cd2 is None:
        raise Exception("Missing CD or CDELT keywords in header")
    wcs = pywcs.WCS(head)

    if units == 'wcs':
        if coordsys=='celestial' and wcs.wcs.lngtyp=='GLON':
            xc,yc = coords.Position((xc,yc),system=coordsys).galactic()
        elif coordsys=='galactic' and wcs.wcs.lngtyp=='RA':
            xc,yc = coords.Position((xc,yc),system=coordsys).j2000()

    if useMontage and CanUseMontage:
        head['CRVAL1'] = xc
        head['CRVAL2'] = yc
        if units == 'pixels':
            head['CRPIX1'] = xw
            head['CRPIX2'] = yw
            head['NAXIS1'] = int(xw*2)
            head['NAXIS2'] = int(yw*2)
        elif units == 'wcs':
            
            cdelt = numpy.sqrt(cd1**2+cd2**2)
            head['CRPIX1'] = xw   / cdelt
            head['CRPIX2'] = yw   / cdelt
            head['NAXIS1'] = int(xw*2 / cdelt)
            head['NAXIS2'] = int(yw*2 / cdelt)

        head.toTxtFile('temp_montage.hdr',clobber=True)
        newfile = montage.wrappers.reproject_hdu(file[0],header='temp_montage.hdr',exact_size=True)
        os.remove('temp_montage.hdr')
    else:

        xx,yy = wcs.wcs_world2pix(xc,yc,0)
        
        if units=='pixels':
            xmin,xmax = numpy.max([0,xx-xw]),numpy.min([head['NAXIS1'],xx+xw])
            ymin,ymax = numpy.max([0,yy-yw]),numpy.min([head['NAXIS2'],yy+yw])
        elif units=='wcs':
            xmin,xmax = numpy.max([0,xx-xw/numpy.abs(cd1)]),numpy.min([head['NAXIS1'],xx+xw/numpy.abs(cd1)])
            ymin,ymax = numpy.max([0,yy-yw/numpy.abs(cd2)]),numpy.min([head['NAXIS2'],yy+yw/numpy.abs(cd2)])
        else:
            raise Exception("Can't use units %s." % units)

        if xmax < 0 or ymax < 0:
            raise ValueError("Max Coordinate is outside of map: %f,%f." % (xmax,ymax))
        if ymin >= head.get('NAXIS2') or xmin >= head.get('NAXIS1'):
            raise ValueError("Min Coordinate is outside of map: %f,%f." % (xmin,ymin))

       
        head['CRPIX1']-=xmin
        head['CRPIX2']-=ymin
        head['NAXIS1']=int(xmax-xmin)
        head['NAXIS2']=int(ymax-ymin)

        if head.get('NAXIS1') == 0 or head.get('NAXIS2') == 0:
            raise ValueError("Map has a 0 dimension: %i,%i." % (head.get('NAXIS1'),head.get('NAXIS2')))

        img = file[0].data[ymin:ymax,xmin:xmax]
        newfile = pyfits.PrimaryHDU(data=img,header=head)
        if verbose: print "Cut image %s with dims %s to %s.  xrange: %f:%f, yrange: %f:%f" % (filename, file[0].data.shape,img.shape,xmin,xmax,ymin,ymax)

    if isinstance(outfile,str):
        newfile.writeto(outfile,clobber=clobber)

    if opened:
        file.close()

    return newfile  
#End of galaxy cutout
