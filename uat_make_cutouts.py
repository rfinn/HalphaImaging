#! /usr/bin/env python 

"""

GOAL:

    to create postage-stamp cutouts of galaxies


PROCEDURE:


    read catalog containing RA, Dec and some measure of galaxy size (e.g. NSA THETA_50)
    keep RA, DEC for galaxies that are on image
    feed results RA, DEC and Size to cutout2d

INPUTS:

    catalog - catalog or subcatalog of NASA Sloan atlas
    halpha filter number - this is the filter number used to refer to the Halpha filter for your group/cluster (e.g. Ha8 would be 8)

USAGE:

    uat_make_cutouts.py catalog image Halpha_filter_number

    To run type CutoutGenerator.py catalog filter_number into the commandline.
    Catalog is the path to the AGC table 
    The filter number is obtainted from the Kitt Peak National observatory webpage.

EXAMPLE:

    python uat_make_cutouts.py full_path_to_catalog input_image halpha_filter_number

REQUIRED MODULES:

NOTES:
    You may need to change the variables Haimage, Rimage and Hawcimage if your image names are different


REFERENCES:

http://docs.astropy.org/en/stable/nddata/utils.html

"""




import numpy as np
from matplotlib import pyplt as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.nddata.utils import Cutout2D
import argparse


parser = argparse.ArgumentParser(description ='Get cutouts for ')
parser.add_argument('--filestring', dest = 'filestring', default = 'hcftr*o00.fits', help = 'string to use to get input files (default = "hcftr*o00.fits")')
parser.add_argument('--s', dest = 's', default = False, action = 'store_true', help = 'Run sextractor to create object catalogs')
parser.add_argument('--c', dest = 'c', default = False, action = 'store_true', help = 'Run scamp')
parser.add_argument('--w', dest = 'w', default = False, action = 'store_true', help = 'Run swarp to create mosaics')
parser.add_argument('--l', dest = 'l', default = False, help = 'List of images to input to swarp')
parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of default config files')
args = parser.parse_args()


catalog=sys.argv[1]
filternumber=sys.argv[2]

pixelscale=.423



#dictionary of Halpha filters
lmin={'4':6573., '8':6606.,'12':6650.,'16':6682.}
lmax={'4':6669., '8':6703.,'12':6747., '16':6779.}

Zmax=(((lmax[filternumber])/6563.)-1)
Zmin=(((lmin[filternumber])/6563.)-1)
print Zmax, Zmin

#Name of images
Haimage=('Hacs_final.fits')
Rimage= ('R_final.fits')
Hawcimage=('Ha_final.fits')
#Hawc is image with continuum 


def get_cutout(image,RA,Dec,size):
    print 'getting cutouts!'




    
    
#Start of Galaxy cutout code 
"""
======
Cutout
======

Generate a cutout image from a .fits file
"""
try:
    import astropy.io.fits as pyfits
    import astropy.wcs as pywcs
except ImportError:
    import pyfits
    import pywcs
import numpy
try:
    import coords
except ImportError:
    pass # maybe should do something smarter here, but I want agpy to install...
try:
    import montage
    import os
    CanUseMontage=True
except ImportError:
    CanUseMontage=False
except Exception:
    CanUseMontage=False

class DimensionError(ValueError):
    pass



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

def makecuts(image,imagefilter):
    hdulist= fits.open(catalog)
    fdata=hdulist[1].data
    print 'Cutting out', image    
    
    voptflag=fdata.VOPT>.1
    vall=fdata.VOPT+~voptflag*fdata.V21
    zall= vall/(3.e5)
    zFlag = (zall > Zmin)&(zall < Zmax)
    
    fdulist = fits.open(image)
    prihdr = fdulist[0].header
    t=fdulist[0].data
    n2,n1=t.shape
    

    w= WCS(image)
    px,py = w.wcs_world2pix(fdata.RA,fdata.DEC,1)
    onimageflag=(px < n1) & (px >0) & (py < n2) & (py > 0)
    
    keepflag=zFlag & onimageflag
    RA=fdata.RA[keepflag]
    DEC=fdata.DEC[keepflag]
    A100=fdata.A100[keepflag]
    AGCNUMBER=fdata.AGCNUMBER[keepflag]
    print 'number of galaxies to keep = ', sum(keepflag)
    
    for i in range(len(RA)):
                
        a=(A100[i]/200.)*60.
        height= ((6*a)/pixelscale)
        width= ((6*a)/pixelscale)
        

        if (a<.01):
            height=120.
            width=120.
            
            
        outimage=(str(AGCNUMBER[i])+'_'+ imagefilter+".fits")
        cutout(image, RA[i], DEC[i],xw=width ,yw=height , outfile=outimage)
        

makecuts(Haimage,'Ha')
makecuts(Rimage,'R')
makecuts(Hawcimage,'Hawc')
