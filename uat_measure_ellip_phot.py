#!/usr/bin/env python

'''

GOAL:

To develop a non-iraf version of ellipse that will:
- fit ellipse to image
- measure flux in concentric elliptical apertures
- write out flux vs semi-major axis

PROCEDURE:
- fit an elliptical aperture at large radius
- keep PA and ellip fixed
- define apertures of smaller and larger radius
- measure R-band flux in apertures
- measure H-alpha flux in the same apertures 
- plot results, normalizing by the flux in the central aperture


REQUIRED MODULES:

- photutils

REFERENCES:


*******
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

from photutils import EllipticalAperture
#from photutils import CircularAperture
from photutils import aperture_photometry


parser = argparse.ArgumentParser(description = 'This code takes an image, and a SExtractor catalog.  ')
parser.add_argument('--imfile', dest = 'imfile', default='A1367-113394-R',help = 'input image, without fits suffix')
parser.add_argument('--imagepath', dest = 'imagepath', default = '', help = 'path to image.  default is current directory')
parser.add_argument('--plot', dest = 'plot', default = False, action = 'store_true', help = 'generate plot of enclosed flux vs radius using matplotlib.  default is False.')
args = parser.parse_args()


if args.plot:
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    from matplotlib.patches import Ellipse

#image_path = '/Users/rfinn/github/H-alpha_cutouts/test/'
image_path = args.imagepath
#im1='A1367-113394-R'
im1 = args.imfile
image = image_path+im1+'.fits'
# read in image
imdat = fits.getdata(image)
#if args.plot:




    
### get ellipse info from sextractor output
SEcat = image_path+im1+'.cat'
cat = fits.getdata(SEcat,2)
'''
cat.X_IMAGE
cat.Y_IMAGE
cat.ELLIPTICITY
cat.THETA_J2000

'''

#find object closest to center

xdim,ydim = imdat.shape
distance = np.sqrt((cat.X_IMAGE - xdim/2.)**2 + (cat.Y_IMAGE - ydim/2.)**2)
objectID = distance == min(distance)

# find R90 - radius enclosing 90% of flux

R90 = cat.FLUX_RADIUS[objectID][0][0]
rmax = 4*R90

position = [(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0])]
theta = np.radians(90.-cat.THETA_J2000[objectID][0])
a = np.arange(2,rmax,1)
b = (1.-cat.ELLIPTICITY[objectID][0])*a

flux = np.zeros(len(a),'f')

for i in range(len(a)):
    ap = EllipticalAperture(position,a[i],b[i],theta)#,ai,bi,theta) for ai,bi in zip(a,b)]
    phot_table = aperture_photometry(imdat, ap)
    flux[i] = phot_table[0][3]


# plot image with outer ellipse

if args.plot:
    plt.figure()
    plt.imshow(imdat,cmap='gray_r')
    plt.colorbar()
    ax = plt.gca()

    ellipse = Ellipse(xy=(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0]), width=a[-1],height=b[-1],edgecolor='r', fc='None', lw=2, angle=theta+90)
    ax.add_patch(ellipse)

# calculate surface brightness in each aperture

area = np.pi*a*b # area of each ellipse

surface_brightness = np.zeros(len(a),'f')

# first aperture is calculated differently

surface_brightness[0] = flux[0]/area[0]

# out apertures need flux from inner aperture subtracted

for i in range(1,len(a)):

    surface_brightness[i] = (flux[i] - flux[i-1])/(area[i]-area[i-1])

# write out photometry
# radius enclosed flux
outfile = open(im1+'_phot.dat','w')#used to be _phot.dat, but changing it to .dat so that it can be read into code for ellipse profiles

outfile.write('# X_IMAGE Y_IMAGE ELLIPTICITY THETA_J2000 \n')
outfile.write('# %.2f %.2f %.2f %.2f \n'%(cat.X_IMAGE[objectID][0],cat.Y_IMAGE[objectID][0],cat.ELLIPTICITY[objectID][0],cat.THETA_J2000[objectID][0]))
outfile.write('# radius enclosed_flux \n')
for i in range(len(a)):
    outfile.write('%.3f %.3f %.3f \n'%(a[i],flux[i],surface_brightness[i]))
outfile.close()

if args.plot:
    plt.figure()
    plt.plot(a,flux,'bo')
    plt.xlabel('semi-major axis (pixels)')
    plt.ylabel('Enclosed flux')
    plt.show()
               


    
