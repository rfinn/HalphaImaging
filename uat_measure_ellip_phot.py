#!/usr/bin/env python

'''

GOAL:

To develop a non-iraf version of ellipse that will:
- fit ellipse to image
- measure flux in concentric elliptical aperture
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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings

from photutils import EllipticalAperture
from photutils import CircularAperture
from photutils import aperture_photometry

image_path = '/Users/rfinn/github/H-alpha_cutouts/test/'
im1='A1367-113394-R'
#A1367-113078.fits

image = image_path+im1+'.fits'
# read in image
imdat = fits.getdata(image)
plt.imshow(imdat)

### get ellipse info from sextractor output
SEcat = image_path+im1+'.cat'
cat = fits.getdata(SEcat,2)
'''
cat.X_IMAGE
cat.Y_IMAGE
cat.ELLIPTICITY
cat.THETA_J2000

'''

'''
# create elliptical apertures

imdata = np.ones((100,100))




positions = [(30., 30.), (40., 40.)]
positions = [(30.,30.)]
apertures = CircularAperture(positions, r=3.)
#apertures1 = EllipticalAperture(positions, 10, 8, np.radians(65.))
phot_table = aperture_photometry(imdat, apertures)
print(phot_table)


print 'got here'

positions = [(30., 30.)]
radii = [3., 4., 5.]
apertures2 = []
for i in range(1):
    apertures2.append(CircularAperture(positions,r=radii[i]))
#apertures2 = [CircularAperture(positions, r=r) for r in radii]
#apertures2 = CircularAperture(positions,
phot_table2 = aperture_photometry(imdata, apertures2[0])
print(phot_table2)

print 'got here!'

'''

position = [(cat.X_IMAGE[0],cat.Y_IMAGE[0])]
theta = np.radians(90.-cat.THETA_J2000)
a = np.arange(2,15,.5)
b = (1.-cat.ELLIPTICITY[0])*a


flux = np.zeros(len(a),'f')

for i in range(len(a)):
    ap = EllipticalAperture(position,a[i],b[i],theta)#,ai,bi,theta) for ai,bi in zip(a,b)]
    phot_table = aperture_photometry(imdat, ap)
    flux[i] = phot_table[0][0]

# calculate surface brightness in each aperture

area = np.pi*a*b # area of each ellipse

surface_brightness = np.zeros(len(a),'f')

# first aperture is calculated differently

surface_brightness[0] = flux[0]/area[0]

# out apertures need flux from inner aperture subtracted

for i in range(1,len(a)):

    surface_brightness[i] = (flux[i] - flux[i-1])/(area[i]-area[i-1])

plt.figure()
plt.plot(a,surface_brightness)
plt.xlabel('semi-major axis (pixels)')
plt.ylabel('Intensity (ADU/pixels^2)')
plt.gca().set_yscale('log')
#plt.gca().set_xscale('log')
plt.savefig(im1+'_sb.png')
#phot_table = aperture_photometry(imdat, apertures)
#print(phot_table['aperture_sum']) 


# fit ellipse
'''
p_init = models.Polynomial2D(degree=2)
fit_p = fitting.LevMarLSQFitter()

x=np.arange(imdat.shape[0])
y=np.arange(imdat.shape[1])
with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
    warnings.simplefilter('ignore')
    p = fit_p(p_init, x,y,imdat)

'''

'''
np.random.seed(0)
y, x = np.mgrid[:imdat.shape[0], :imdat.shape[1]]
#z = 2. * x ** 2 - 0.5 * x ** 2 + 1.5 * x * y - 1.
#z += np.random.normal(0., 0.1, z.shape) * 50000.
z = imdat
# Fit the data using astropy.modeling
x0 = imdat.shape[1]/2.
y0 = imdat.shape[0]/2.
a = imdat.shape[0]/5.
b = .6*a
theta=Angle(0,'deg')
p_init = models.Ellipse2D(amplitude=np.max(imdat),x_0 = x0,y_0 = y0, a=a, b=b, theta=theta.radian)
fit_p = fitting.LevMarLSQFitter()

#imdat=50000*imdat
with warnings.catch_warnings():
    # Ignore model linearity warning from the fitter
    warnings.simplefilter('ignore')
    #p = fit_p(p_init, x, y, z)
    p = fit_p(p_init, x, y, imdat)

# Plot the data with the best-fit model
plt.figure(figsize=(8, 2.5))
plt.subplot(1, 3, 1)
plt.imshow(imdat, origin='lower', interpolation='nearest', vmin=-1, vmax=12)
plt.title("Data")
plt.subplot(1, 3, 2)
plt.imshow(p(x, y), origin='lower', interpolation='nearest', vmin=-1, vmax=12)
plt.title("Model")
plt.subplot(1, 3, 3)
plt.imshow(imdat - p(x, y), origin='lower', interpolation='nearest', vmin=-1,vmax=12)
plt.title("Residual")
# display image and ellipse
'''
