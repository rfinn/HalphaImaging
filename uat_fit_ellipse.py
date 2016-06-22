#!/usr/bin/env python

'''


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
