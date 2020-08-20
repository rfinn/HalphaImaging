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

example of fitting an elliptical aperture at 

http://docs.astropy.org/en/stable/api/astropy.modeling.functional_models.Ellipse2D.html#astropy.modeling.functional_models.Ellipse2D

import numpy as np
from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
x0, y0 = 25, 25
a, b = 20, 10
theta = Angle(30, 'deg')
e = Ellipse2D(amplitude=100., x_0=x0, y_0=y0, a=a, b=b,
              theta=theta.radian)
y, x = np.mgrid[0:50, 0:50]
fig, ax = plt.subplots(1, 1)
ax.imshow(e(x, y), origin='lower', interpolation='none', cmap='Greys_r')
e2 = mpatches.Ellipse((x0, y0), 2*a, 2*b, theta.degree, edgecolor='red',
                      facecolor='none')
ax.add_patch(e2)
plt.show()


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
from astropy.modeling import models, fitting
from astropy.modeling.models import Ellipse2D
from astropy.coordinates import Angle
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings

image_path = '/Users/rfinn/github/H-alpha_cutouts/R-band-cutouts/'
im1='A1367-113394-R.fits'
#A1367-113078.fits

image = image_path+im1
# read in image
imdat = fits.getdata(image)
plt.imshow(imdat)



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
