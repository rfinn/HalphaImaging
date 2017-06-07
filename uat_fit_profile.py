#!/usr/bin/env python




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


'''

cat1 =  'A1367-113068-R_phot.dat'
cat2 =  'A1367-113068-CS_phot.dat'
R = loadtxt(cat1)
R
CS = loadtxt(cat2)
figure()
clf()
plot(CS[:,0],CS[:,2]/CS[:,2][0],'b.')
plot(R[:,0],R[:,2]/R[:,2][0],'r.')
gca().set_yscale('log')


'''


# read in intensity versus radius


# calculate surface


plt.figure()
plt.plot(a,surface_brightness)
plt.xlabel('semi-major axis (pixels)')
plt.ylabel('Intensity (ADU/pixels^2)')
plt.gca().set_yscale('log')
#plt.gca().set_xscale('log')
plt.savefig(im1+'_sb.png')
#phot_table = aperture_photometry(imdat, apertures)
#print(phot_table['aperture_sum']) 




# fit exponential profile

