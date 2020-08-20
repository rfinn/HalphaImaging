#!/usr/bin/env python

'''

GOAL:

To develop a non-iraf version of ellipse that will:
- fit ellipse to image
- measure flux in concentric elliptical apertures
- write out flux vs semi-major axis

PROCEDURE:
- get PA and ellip from SExtractor
- keep PA and ellip fixed
- define apertures of smaller and larger radius
- measure R-band flux in apertures
- measure H-alpha flux in the same apertures 
- plot results, normalizing by the flux in the central aperture

EXAMPLE:

in ipython:

%run ~/github/HalphaImaging/uat_measure_ellip_phot.py --imfile pointing-4-88353-R --mask pointing-4-88353-mask.fits --plot
%run ~/github/HalphaImaging/uat_measure_ellip_phot.py --imfile pointing-4-88353-CS --mask pointing-4-88353-mask.fits --plot


REQUIRED MODULES:

- photutils

REFERENCES:


*******
NOTES ON HOW TO USE photutils

Measure aperture flux

http://photutils.readthedocs.io/en/latest/photutils/aperture.html


>>> from photutils import EllipticalAperture
>>> from photutils import aperture_photometry
>>> a = 5.
>>> b = 3.
>>> theta = np.pi / 4.
>>> apertures = EllipticalAperture(positions, a, b, theta)
>>> phot_table = aperture_photometry(data, apertures)
>>> print(phot_table['aperture_sum'])    


'''

from astropy.io import fits
import numpy as np
#from astropy.modeling import models, fitting
#from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
import warnings
import argparse
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
from photutils import EllipticalAperture
#from photutils import CircularAperture
from photutils import aperture_photometry



parser = argparse.ArgumentParser(description = 'This code takes an image, and a SExtractor catalog.  ')
parser.add_argument('--pointing',dest = 'pointing', default = None, help = 'Cluster prefix of image for continuum subtraction.  Use this if you are subtracting continuum from a cutout rather than a mosaic.')
parser.add_argument('--id',dest = 'id', default = None, help = 'NSAID of image for continuum subtraction.  Use this if you are subtracting continuum from a cutout rather than a mosaic.')
#parser.add_argument('--r', dest = 'r', default='A1367-113394-R.fits',help = 'r-band image')
#parser.add_argument('--ha', dest = 'ha', default='A1367-113394-R.fits',help = 'ha image')
parser.add_argument('--mask', dest = 'mask', default=None,help = 'mask to use when measuring photometry')
parser.add_argument('--imagepath', dest = 'imagepath', default = '', help = 'path to image.  default is current directory')
parser.add_argument('--rmax', dest = 'rmax', default = 2., help = 'maximum radius to measure photometry out to.  This is a multiple of what sextractor measures as R90 for r-band image.  Default is 2xR90.')
parser.add_argument('--naper', dest = 'naper', default = 20, help = 'Number of apertures to use.  Default is 20.')
parser.add_argument('--plot', dest = 'plot', default = False, action = 'store_true', help = 'generate plot of enclosed flux vs radius using matplotlib.  default is False.')
parser.add_argument('--verbose', dest = 'verbose', default = False, action = 'store_true', help = 'print statements to monitor progress of code')
args = parser.parse_args()



if args.plot:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import Ellipse

#image_path = '/Users/rfinn/github/H-alpha_cutouts/test/'
image_path = args.imagepath

if args.id:
    name_r = args.pointing+'-'+args.id+'-r.fits'
    name_ha = args.pointing+'-'+args.id+'-CS.fits'
    id = args.id


rimage = image_path+name_r
# read in image
imdat_r, header_r = fits.getdata(rimage, header=True)

haimage = image_path+name_ha
# read in image
imdat_ha, header_ha = fits.getdata(haimage, header=True)

print('reading sextractor catalog\n')    
### get ellipse info from sextractor output
SEcat = image_path+name_r.split('.fits')[0]+'.cat'
cat = fits.getdata(SEcat,2)
'''
cat.X_IMAGE
cat.Y_IMAGE
cat.ELLIPTICITY
cat.THETA_J2000

'''

#find object closest to center

xdim,ydim = imdat_r.shape
distance = np.sqrt((cat.X_IMAGE - xdim/2.)**2 + (cat.Y_IMAGE - ydim/2.)**2)
objectID = distance == min(distance)

# find R90 - radius enclosing 90% of flux

R90 = cat.FLUX_RADIUS[objectID][0][0]
rmax = float(args.rmax)*R90

print('max radius for measuring photometry is '+str(rmax))
position = [(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0])]

## create an array of central positions
#positions = np.zeros([int(args.naper),2],'f')
#positions[:,0] = np.ones(int(args.naper))*cat.X_IMAGE[objectID][0]
#positions[:,1] = np.ones(int(args.naper))*cat.Y_IMAGE[objectID][0]

# you can't define multiple concentric apertures
# because a, b, and theta can't be arrays :(

# define ellipse shape
theta = np.radians(90.-cat.THETA_J2000[objectID][0])
a = np.linspace(2,rmax,int(args.naper))
b = (1.-cat.ELLIPTICITY[objectID][0])*a

print('ellipticity = ',cat.ELLIPTICITY[objectID][0])
print('b/a = ',1.-cat.ELLIPTICITY[objectID][0])

flux_r = np.zeros(len(a),'f')
flux_ha = np.zeros(len(a),'f')
if args.mask:
    maskdat = fits.getdata(args.mask)
    maskdat = np.array(maskdat,'bool')

for i in range(len(a)):
    if args.verbose:
        print('defining elliptical aperture '+str(i)+' ap size = ',str(a[i]))
    ap = EllipticalAperture(position,a[i],b[i],theta)#,ai,bi,theta) for ai,bi in zip(a,b)]

    if args.mask:
        phot_table_r = aperture_photometry(imdat_r, ap, mask=maskdat)
        phot_table_ha = aperture_photometry(imdat_ha, ap, mask=maskdat)
    else:
        phot_table_r = aperture_photometry(imdat_r, ap)
        phot_table_ha = aperture_photometry(imdat_ha, ap)
    if args.verbose:
        print('measuring flux in aperture '+str(i)+' ap size = ',str(a[i]))
    flux_r[i] = phot_table_r['aperture_sum'][0]
    flux_ha[i] = phot_table_ha['aperture_sum'][0]

# plot image with outer ellipse

if args.plot:
    if args.verbose:
        print('plotting results \n')
    plt.figure(figsize=(10,4))
    plt.subplot(1,2,1)
    
    vmin,vmax=scoreatpercentile(imdat_r,[.5,99.5])
    plt.imshow(imdat_r,cmap='gray_r',vmin=vmin,vmax=vmax,origin='lower')
    plt.colorbar()
    ax = plt.gca()

    ellipse = Ellipse(xy=(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0]), width=2*a[-1],height=2*b[-1],edgecolor='r', fc='None', lw=2, angle=-1*np.degrees(theta))
    ax.add_patch(ellipse)

    plt.subplot(1,2,2)
    vmin,vmax=scoreatpercentile(imdat_ha,[.5,99.5])
    plt.imshow(imdat_ha,cmap='gray_r',vmin=vmin,vmax=vmax,origin='lower')
    plt.colorbar()
    ax = plt.gca()

    ellipse = Ellipse(xy=(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0]), width=2*a[-1],height=2*b[-1],edgecolor='r', fc='None', lw=2, angle=-1*np.degrees(theta))
    ax.add_patch(ellipse)
    
    plt.savefig(name_r.split('.fits')[0]+'-snapshot.png')
# calculate surface brightness in each aperture

area = np.pi*a*b # area of each ellipse

surface_brightness_r = np.zeros(len(a),'f')
surface_brightness_ha = np.zeros(len(a),'f')

# first aperture is calculated differently

surface_brightness_r[0] = flux_r[0]/area[0]
surface_brightness_ha[0] = flux_ha[0]/area[0]

# out apertures need flux from inner aperture subtracted

for i in range(1,len(a)):

    surface_brightness_r[i] = (flux_r[i] - flux_r[i-1])/(area[i]-area[i-1])
    surface_brightness_ha[i] = (flux_ha[i] - flux_ha[i-1])/(area[i]-area[i-1])

# write out photometry for r-band
# radius enclosed flux
outfile = open(name_r.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles

outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
outfile.write('# %.2f %.2f %.2f %.2f \n'%(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0],cat.ELLIPTICITY[objectID][0],cat.THETA_J2000[objectID][0]))
outfile.write('# radius enclosed_flux \n')
for i in range(len(a)):
    outfile.write('%.2f %.3e %.3e \n'%(a[i],flux_r[i],surface_brightness_r[i]))
outfile.close()

# write out photometry for h-alpha
# radius enclosed flux
outfile = open(name_ha.split('.fits')[0]+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles

outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
outfile.write('# %.2f %.2f %.2f %.2f \n'%(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0],cat.ELLIPTICITY[objectID][0],cat.THETA_J2000[objectID][0]))
outfile.write('# radius enclosed_flux \n')
for i in range(len(a)):
    outfile.write('%.2f %.3e %.3e \n'%(a[i],flux_ha[i],surface_brightness_ha[i]))
outfile.close()

if args.plot:
    plt.figure(figsize=(10,4))
    plt.subplots_adjust(wspace=.3)
    plt.subplot(1,2,1)
    plt.plot(a,flux_r,'bo')
    plt.xlabel('semi-major axis (pixels)')
    plt.ylabel('Enclosed flux')
    plt.title(name_r)
    plt.subplot(1,2,2)
    plt.plot(a,flux_ha,'bo')
    plt.xlabel('semi-major axis (pixels)')
    plt.ylabel('Enclosed flux')
    plt.title(name_ha)
    plt.show()
               


    
