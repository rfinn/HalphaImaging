#!/usr/bin/env python

'''
USAGE:

from within ipython:
%run ~/github/HalphaImaging/getzp.py --image pointing031-r.coadd.fits --instrument i --filter r

The y intercept is -1*ZP

To print the value in ipython, type:

-1*zp.bestc[1]


To run on mosaic data, set the flatten and spline flags:

%run ~/github/HalphaImaging/getzp.py --image pointing031-r.coadd.fits --instrument m --filter r --flatten 1 --spline

UPDATES:
* 2026-Jan: updating color transformations
* implemented scipy.optimize.curve_fit in getzp.py to 
    * keep slope fixed at 1
    * get an estimate of error in ZP (sqrt(covariance))
* program now prints ZP and error at the end

NOTES:

2019-06-12
* making sure saturated stars are ignored
- coadd produced by swarp is in ADU/exptime
- added argument nexptime that allows user to toggle between images in ADU vs ADU/s.  If image is in ADU/s, then I grab the exptime from the image header and change SATUR_LEVEL to 40000./exptime


Apertures:
- by default, we are using aperture magnitudes, aperture 5 is default, which is
REFERENCES:

Pan-STARRS
https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/

https://panstarrs.stsci.edu/


GAIA
https://gea.esac.esa.int/archive/documentation/GDR1/Data_processing/chap_cu5phot/sec_phot_calibr.html

https://www.cosmos.esa.int/web/gaia/dr2-known-issues


OLD SDSS QUERY
from astroquery.sdss import SDSS


query = 'SELECT TOP 10 ra, dec, u,g,r,i,z, flags_r FROM Star WHERE (clean = 1) AND ra BETWEEN 180 and 181 AND dec BETWEEN -0.5 and 0.5 AND ((flags_r & 0x10000000) != 0) AND ((flags_r & 0x8100000c00a4) = 0) AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.2)) AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)'

t = SDSS.query_sql(query, data_release=14)
'''


import argparse
import os
import numpy as np
import sys
import glob
import time
start_time = time.time()

from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib
#matplotlib.use("Tkagg")
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.stats import mad_std
from astropy.table import Table, hstack, Column
from astropy.stats import sigma_clip
from astropy.modeling import models, fitting

from scipy.optimize import curve_fit
from scipy import interpolate
try:
    from scipy.stats import median_abs_deviation as MAD2
except:
    from scipy.stats import median_absolute_deviation as MAD2
from astroquery.vizier import Vizier
try:
    from photutils import Background2D, MedianBackground # photutils v 1.3.0
except ImportError: # photutils v 2.2
    from photutils.background import Background2D, MedianBackground
import itertools

#################################################
## GLOBAL VARIABLES
#################################################
FIT_SPLINE = True

# function for fitting ZP equation
# this one forces slope = 1
zpfunc = lambda x, zp: x + zp

# this function allows the slope to vary
zpfuncwithslope = lambda x, m, zp: m*x + zp

pixelscale = {'HDI':0.43, 'INT':0.333, 'BOK':0.45252,'MOS':0.425 } 

def panstarrs_query(ra_deg, dec_deg, rad_deg, maxmag=19,
                    maxsources=10000):
    """
    FOUND THIS ON THE WEB 
    https://michaelmommert.wordpress.com/2017/02/13/accessing-the-gaia-and-pan-starrs-catalogs-using-python/

    
    Query PanSTARRS @ VizieR using astroquery.vizier
    :param ra_deg: RA in degrees
    :param dec_deg: Declination in degrees
    :param rad_deg: field radius in degrees
    :param maxmag: upper limit G magnitude (optional)
    :param maxsources: maximum number of sources
    :return: astropy.table object

    https://vizier.cds.unistra.fr/viz-bin/VizieR-2
    """
    pan_columns =['objID', 'RAJ2000', 'DEJ2000','e_RAJ2000', 'e_DEJ2000', 'f_objID', 'Qual','gmag', 'e_gmag','rmag', 'e_rmag','imag', 'e_imag','zmag', 'e_zmag','ymag', 'e_ymag']
    #print(pan_columns)
    
    vquery = Vizier(columns=pan_columns,column_filters={"gmag":("<%f" % maxmag)},row_limit=maxsources)
    print('HEY!!! in panstarrs_query, ra,dec,rad = ',ra_deg,dec_deg,rad_deg)
    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg, unit=(u.deg, u.deg),frame='icrs')
    t = vquery.query_region(field, width=("%fd" % rad_deg), catalog="II/349/ps1")
    print("in panstarrs_query...")
    print(t)
    return t[0]

def fitline_sigma_clipping(x, y, yerr=None, nsigma=3,niter=3):
    """ 
    fit a line with slope fixed at 1  

    PARAMS:
    * x
    * y
    * yerr OPTIONAL

    RETURNS:
    * fit function
    * mask = with rejected values=True
    * uncertainty = np array with uncertainty in slope

    """
    # updating to get initial estimate of intercept from polyfit
    # then feed this into Linear1D
    # hoping this helps the more difficult cases converge
    c = np.polyfit(x,y,1)
    intercept = c[1]
    print(f"initial fit from polyfit: intercept={intercept:.2f}")
    # 2. Initialize the polynomial model and fitter
    poly_init = models.Linear1D(intercept=intercept, fixed={'slope': True},slope=1) # Initialize a 2nd degree polynomial model

    fitter = fitting.LinearLSQFitter() # Standard least squares fitter
    
    # 3. Initialize the outlier removal fitter
    # This wraps the standard fitter with sigma clipping logic
    or_fit = fitting.FittingWithOutlierRemoval(
        fitter,
        sigma_clip,
        niter=niter,        # Number of iterations for sigma clipping
        sigma=nsigma       # Number of standard deviations to clip beyond
    )
    
    # 4. Fit the data
    # The fit returns the fitted model and a mask indicating clipped points
    if yerr is not None:
        fitted_poly, mask = or_fit(poly_init, x, y, weights=1/yerr)
    else:
        fitted_poly, mask = or_fit(poly_init, x, y)
    
    return fitted_poly, mask, or_fit.fitter.fit_info['singular_values']

def fitspline2d(x,y,z,nx,ny,order=3,s=1000):
    """
    Fit a 2d spline to a surface given a set of (x,y,z) coordinates

    I am using this to fit a surface to the residuals of measures vs panstarrs flux,
    so this is basically making an flatfield image.

    INPUT:
    x,y = positions of points
    z = in this case, z is the ratio of f_obs/f_panstarrs
    nx = dimension of image in x direction
    ny = dimension of image in x direction
    norder = 3; order of the spline
    s = 1000; smoothing factor for spline

    RETURN:
    zim = 2d image (ny,nx) that represents the spline fit; you can normalize input image by this
          to obtain "flattened" image

    REFERENCE:
    https://docs.scipy.org/doc/scipy/tutorial/interpolate/smoothing_splines.html#d-smoothing-splines

    """
    print()
    print("fitting 2d spline.  go refresh your coffee...\n")
    xnew_edges, ynew_edges = np.mgrid[0:nx+1:complex(nx+1), 0:ny+1:complex(ny+1)]
    xnew = xnew_edges[:-1, :-1] + np.diff(xnew_edges[:2, 0])[0] / 2.
    ynew = ynew_edges[:-1, :-1] + np.diff(ynew_edges[0, :2])[0] / 2.
    tck = interpolate.bisplrep(x, y, z, s=s,kx=order,ky=order)
    znew = interpolate.bisplev(xnew[:,0], ynew[0,:], tck)
    zim = np.transpose(znew)
    # testing to see if transpose is not necessary
    #zim = znew
    print("returning to your regular program...")
    return zim
    
def polyfit2d(x, y, z, order=3):
    # from  https://stackoverflow.com/questions/7997152/python-3d-polynomial-surface-fit-order-dependent
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, _, _, _ = np.linalg.lstsq(G, z)
    return m

def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    #z = np.zeros_like(x)
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.zeros(x.shape, dtype=float)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

def fit_circle_func():
    """
    #fit function centered on center of FOV
    never developed this...
    """
    pass

def get_filebasename(fname):
    # determine telescope and filter from filename
    t = fname.split('-')
    #print(t)
    # check to see if declination is negative
    # if it is, then index of instrument will be off by one
    try:
        if fname[10] == '+':
            basename = '-'.join(t[0:5])
        else:
            basename = '-'.join(t[0:6])
    except IndexError:
        basename = fname
    return basename
        

class args():
    """ for replicating argparse input to getzp without using argparse"""
    def __init__(self,image,instrument,filter,nexptime=False):
        self.catdir = 'SEcats_getzp'
        self.image = image
        self.instrument = instrument
        self.filter = filter
        self.normbyexptime = nexptime
        self.verbose = False
        self.nsigma = 3.
        
        self.d = os.getenv("HOME")+'/github/HalphaImaging/astromatic/'
        self.catalog = None
        self.fwhm = None
        self.useri = False
        self.mag = 0
        self.naper = 5
        self.fit = False
        self.flatten = 0
        self.norder = 2
        self.fixbok = False
        
class getzp():
    def __init__(self, args):
        self.catdir = 'SEcats_getzp'
        self.image = args.image

        header = fits.getheader(self.image)
        self.header = header
        #self.plotprefix = self.image.split('coadd')[0].replace('.','-').replace('-noback',"")
        #try:
        #    self.plotprefix = header['OBJECT']
        #except KeyError:
        #    self.plotprefix = self.image.split('.fits')[0]
        self.plotprefix = self.image.split('.fits')[0]+'-'
            
        # create plot directory if it doesn't already exist
        if not os.path.exists('plots'):
            os.mkdir('plots')
        self.verbose = args.verbose
        self.astrodir = args.d
        self.instrument = args.instrument
        if self.instrument == 'h':
            self.pixelscale = pixelscale['HDI']
        elif self.instrument == 'i':
            self.pixelscale = pixelscale['INT']
        elif self.instrument == 'b':
            self.pixelscale = pixelscale['BOK']            
        elif self.instrument == 'm':
            self.pixelscale = pixelscale['MOS']            
        self.filter = args.filter

        self.MINMAG = float(args.minmag)
        #if self.instrument == 'i':
        #    self.MINMAG = 15.
        #if self.instrument == 'b':
        #    self.MINMAG = 15.5
        # if image is in ADU rather than ADU/s, divide image by exptime before running sextractor
        if args.normbyexptime:
            im, header = fits.getdata(self.image,header=True)
            exptime = header['EXPTIME']
            norm_im = im/float(exptime)
            header.set('ORIGEXPT',value=exptime,comment='ORIGINAL EXPTIME BEFORE NORM')
            header.set('EXPTIME',value=1,comment="getzp set exptime=1")
            fits.writeto('n'+self.image, norm_im, header, overwrite=True)
            self.image = 'n'+self.image
            self.plotprefix = 'n'+self.plotprefix

        if self.verbose:
            print('output image = ', self.image)
        self.nsigma = float(args.nsigma)
        self.useri = args.useri
        self.naper = int(args.naper)
        self.mag = int(args.mag)
        self.flatten = int(args.flatten)
        self.spline = args.spline
        self.spline_order = int(args.spline_order)
        self.spline_smooth = int(args.spline_smooth)
        self.norder = int(args.norder)
        self.fwhm = args.fwhm
        self.fixbok = args.fixbok
        self.fixint = args.fixint
        self.minmag = args.minmag
        self.args = args
        global v1, v2
        if 'ha' in self.filter:
            v1 = .9
            v2 = 1.1
            #v1 = .95
            #v2 = 1.05 #  keep halph range to match broad band filters
        else:
            v1=.95
            v2=1.05


    def getzp(self):
        plt.close('all')
        if self.verbose:
            print('')
            print('STATUS: running se')
        secat = os.path.join(self.catdir,self.image.split('.fits')[0]+'.cat')
        if os.path.exists(secat):
            print("found SE cat!!!")
            self.read_se_cat()
        else:
            self.runse(useprevious=False)
        if self.verbose:
            print('')        
            print('STATUS: getting panstarrs')
        self.get_panstarrs()
        self.plot_se_pan_positions()
        if self.verbose:
            print('')        
            print('STATUS: matching se cat to panstarrs')       
        self.match_coords()
        if self.verbose:
            print('')        
            print('STATUS: fitting zeropoint')        
        self.fitzp(plotall=True)
        if self.verbose:
            print('')        
            print('STATUS: udating header')
        if self.instrument == 'b':
            if self.fixbok:  
                print("checking 90prime ccds...",self.fixbok)
                self.check_90prime_ccds()
            self.runse(useprevious=False)
            self.match_coords()
            self.fitzp(plotall=True)
        elif self.instrument == 'i':
            if self.fixint:
                print("checking INT ccds...", self.fixint)
                self.check_INT_ccds(min_stars=10,max_scale_dev=0.1)
            self.runse(useprevious=False)
            self.match_coords()
            self.fitzp(plotall=True)
            
        self.update_header()

        if self.flatten > 0:
            self.fit_residual_surface(norder=self.norder,suffix=None)
            self.renorm_wfc()
            self.rerun_zp_fit()

            if self.flatten == 2:
                print("running a second round of flattening")
                self.fit_residual_surface(norder=self.norder)
                self.renorm_wfc()
                # this creates 'ff'+imagename
                self.rerun_zp_fit()
    
    def getzp_wfc(self):
        self.getzp()
        #self.fit_residual_surface(norder=2)
        # this creates 'f'+imagename
        if self.flatten > 0:
            self.renorm_wfc()
            self.rerun_zp_fit()

            if self.flatten == 2:
                print("running an additional round of flattening for halpha")
                self.fit_residual_surface(norder=self.norder)
                self.renorm_wfc()
                # this creates 'ff'+imagename
                self.rerun_zp_fit()
        # clean up
        # leaving the extra images for now b/c not sure if flattening is need for 2022 INT data
        # TODO - reapply clean up after situation is resolved with flattening 2022 data
        #if self.image.find('f') > -1:
        #    rootname = self.image.strip('f')
        #    if rootname.startswith('nWFC'):
        #        os.remove(rootname)
        #    if self.image.find('ff') > -1:
        #        os.remove('f'+rootname)
                
        plt.figure()
        plt.hist(self.zim,bins=np.linspace(.9,1.1,40))
        plt.axvline(x=1,color='k')
        plt.xlabel('Flux/Flux_predicted of fitted sources')
        
    def runse(self,useprevious=True):
        """
        Run Source Extractor on image to measure magnitudes
        """

        # check for se catalog

        

        t = self.image.split('.fits')
        froot = t[0]
        # check for se catalog
        secat = f"{self.catdir}/{froot}.cat"
        if (os.path.exists(secat) & useprevious):
            print("Found SE catalog: ",secat)
            print("Not rerunning source extractor...")
            return

        os.system('cp ' +self.astrodir + '/default.* .')        
        if self.instrument == 'h':
            defaultcat = 'default.sex.HDI'
        elif self.instrument == 'i':
            defaultcat = 'default.sex.INT.getzp'
            self.keepsection=[1000,5000,0,4000]
        elif self.instrument == 'm':
            defaultcat = 'default.sex.HDI'
        elif self.instrument == 'b':
            print("hey Rose - ")
            print("using default.sex.BOK!!!")
            print()
            defaultcat = 'default.sex.BOK'#.getzp'
        header = fits.getheader(self.image)
        try:
            expt = header['EXPTIME']
        except KeyError:
            expt = 1.
        ADUlimit = 40000.
        if self.instrument == 'm':
            ADUlimit = 40000.
        catdir = self.catdir
        if not os.path.exists(catdir):
            os.mkdir(catdir)
            
        if self.fwhm is None:
            
            t = f'sex {self.image} -c {defaultcat} -CATALOG_NAME {self.catdir}/{froot}.cat -BACK_SIZE 128 -BACKPHOTO_TYPE LOCAL -BACK_TYPE AUTO -MAG_ZEROPOINT 0 -SATUR_LEVEL {ADUlimit:.1f}'
            if self.instrument == 'b':
                weight_image = self.image.replace('.fits','.weight.fits')
                if weight_image.startswith('f'):
                    weight_image = weight_image[1:]
                if os.path.exists(weight_image):
                    t += f" -WEIGHT_IMAGE {weight_image} -WEIGHT_TYPE MAP_WEIGHT -RESCALE_WEIGHTS  N"
                else:
                    print(f"WARNING: no weight image {weight_image}")
            if self.verbose:
                print('running SE first time to get estimate of FWHM')
                print(t)
            os.system(t)

            # clean up SE files
            # skipping for now in case the following command accidentally deletes user files
            # os.system('rm default.* .')


        ###################################
        # Read in Source Extractor catalog
        ###################################
        if self.verbose:
            print('reading in SE catalog from first pass')
        secat_filename = f"{catdir}/{froot}.cat"
        self.secat = fits.getdata(secat_filename,2)
        self.secat0 = self.secat
        # get median fwhm of image
        # for some images, this comes back as zero, and I don't know why
        fwhm = np.median(self.secat['FWHM_IMAGE'])*self.pixelscale
            
            
        t = f'sex {self.image} -c {defaultcat} -CATALOG_NAME {catdir}/{froot}.cat -BACK_SIZE 128 -BACKPHOTO_TYPE LOCAL -BACK_TYPE AUTO -MAG_ZEROPOINT 0 -SATUR_LEVEL {ADUlimit:.1f} -SEEING_FWHM {fwhm:.2f}'
        if self.instrument == 'b':
            weight_image = self.image.replace('.fits','.weight.fits')
            if os.path.exists(weight_image):
                t += f" -WEIGHT_IMAGE {weight_image} -WEIGHT_TYPE MAP_WEIGHT -RESCALE_WEIGHTS  N"
            else:
                print(f"WARNING: no weight image {weight_image}")
        
        if float(fwhm) == 0:
            print('WARNING: measured FWHM is zero!')
            if self.verbose:
                print('running SE again with new FWHM to get better estimate of CLASS_STAR')
        else:
            if self.verbose:
                print(t)
                print('running SE w/user input for FWHM to get better estimate of CLASS_STAR')            
        #############################################################
        # rerun Source Extractor catalog with updated SEEING_FWHM
        #############################################################

        #print(t)
        os.system(t)
        self.read_se_cat()
    def read_se_cat(self):
        

        ###################################
        # Read in Source Extractor catalog
        ###################################
        t = self.image.split('.fits')
        froot = t[0]

        # adding subdir for catalogs
        secat_filename = f"{self.catdir}/{froot}.cat"
        self.secat = fits.getdata(secat_filename,2)

        
        ###################################
        # get max/min RA and DEC for the image
        ###################################

        minRA = min(self.secat['ALPHA_J2000'])
        maxRA = max(self.secat['ALPHA_J2000'])
        minDEC = min(self.secat['DELTA_J2000'])
        maxDEC = max(self.secat['DELTA_J2000'])
        self.centerRA = 0.5*(minRA + maxRA)
        self.centerDEC = 0.5*(minDEC + maxDEC)
        #radius = np.sqrt((maxRA- centerRA)**2 + (maxDEC - centerDEC)**2)
        #print(radius)

        # NOTE - radius is with width of the rectangular
        # search area.  This is not a circular search, as I orginally thought.
        self.radius = np.sqrt(2)*max((maxRA - minRA), (maxDEC - minDEC))/2
        # fixing to be compatible with box search geometry
        self.width = max((maxRA - minRA), (maxDEC - minDEC))
    def get_panstarrs(self):

        ###################################
        # get Pan-STARRS catalog over the same region
        ###################################
        # this naming convention uses the same panstarrs catalog
        # for both the r and halpha images - saves time!
        ptab_name = get_filebasename(self.image).replace('.fits','_pan_tab.csv') # don't see how this works - there is no fits in the basename
        ptab_name = get_filebasename(self.image)+'_pan_tab.csv'
        print(f"self.image={self.image}, panstarrs table = {ptab_name}")
        if os.path.exists(ptab_name):
            print('panstarrs table already downloaded')
            self.pan = Table.read(ptab_name)
            
        else:

            print()
            print("no previous panstarrs catalog found :(")
            self.pan = panstarrs_query(self.centerRA, self.centerDEC, self.width)
            ptab = Table(self.pan)
            ptab.write(ptab_name,format='csv',overwrite=True)
        self.color_correct_panstarrs()
    def match_coords(self):
        ###################################
        # match Pan-STARRS1 data to Source Extractor sources
        # remove any objects that are saturated or non-linear in our r-band image
        ###################################


        secoords = SkyCoord(ra=self.secat['ALPHA_J2000'],dec=self.secat['DELTA_J2000'],unit=(u.degree,u.degree),frame='icrs')
        #print(self.pan)
        pancoords = SkyCoord(ra=self.pan['RAJ2000'],dec=self.pan['DEJ2000'],unit=(u.degree,u.degree),frame='icrs')

        index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)

        # only keep matches with matched RA and Dec w/in 5 arcsec
        self.matchflag = dist2d.degree < 5./3600
        self.matchedarray1=np.zeros(len(pancoords),dtype=self.secat.dtype)
        self.matchedarray1[self.matchflag] = self.secat[index[self.matchflag]]

        ###################################################################
        # remove any objects that are saturated, have FLAGS set, galaxies,
        # must have 15 < r < 16.8 for panstarrs mag, based on fitting results
        ###################################################################
        if self.verbose:
            print(f'\t matched {np.sum(self.matchflag)} objects')

        # TODONE - add color restriction to panstarrs: 0 < g-r < 1
        ps_gr = self.pan['gmag'] - self.pan['rmag']
        colorflag = (ps_gr > 0) & (ps_gr < 1)
        
        #colorflag = (ps_gr > -0.2) & (ps_gr < 1.2)

        # testing what happens with removal of color cut
        # doesn't affect anything
        #colorflag = np.ones(len(ps_gr),'bool')
        self.fitflag = colorflag & self.matchflag  & (self.pan['rmag'] > self.MINMAG) & (self.matchedarray1['FLAGS'] <  1) & (self.pan['Qual'] < 64)  & (self.pan['rmag'] < 19) #& (self.matchedarray1['CLASS_STAR'] > 0.95) #& (self.matchedarray1['MAG_AUTO'] > -11.)
        
        if self.verbose:
            print(f'\t number that pass fit {np.sum(self.fitflag)}')
        # for WFC on INT, restrict area to central region
        # to avoid top chip and vignetted regions
        #
        # skipping for now in place of fitting the residuals
        # and re-normalizing the image as per suggestion of mischa
        #
        #if self.instrument == 'i':
        #    self.goodarea_flag = (self.matchedarray1['X_IMAGE'] > self.keepsection[0]) & #\
        #        (self.matchedarray1['X_IMAGE'] < self.keepsection[1]) & \
        #        (self.matchedarray1['Y_IMAGE'] > self.keepsection[2]) & \
        #        (self.matchedarray1['Y_IMAGE'] < self.keepsection[3])
        #    self.fitflag = self.fitflag & self.goodarea_flag
    def color_correct_panstarrs(self):
        """
        from values that I calculated from filter traces and Manga MaStar normal spectral

        - below are coefficients for second order polynomial fits
        - first value is coefficient of (PS_g-PS_r)^2

        - the values are calculated in the notebook github/filter_transformations/notebooks/filtertrans-dev.ipynb

        """

        # fit using stars with 0 < g-r < 1
        filter_trans_dict = {
            'BOK90prime-r': [0.0092,-0.1076,0.0202,],\
            'BOK90prime-ha4': [0.0336,-0.2664,0.0507,],\
            'MOS-r': [0.0026,-0.0370,0.0055,],\
            'MOS-ha4': [0.0336,-0.2664,0.0507,],\
            'MOS-ha8': [0.0151,-0.2452,0.0311,],\
            'MOS-ha12': [0.0149,-0.2617,0.0404,],\
            'MOS-R': [0.0181,-0.1657,0.0056,],\
            'MOS-ha16': [0.0171,-0.2799,0.0512,],\
            'HDI-r': [0.0013,-0.0095,0.0000,],\
            'HDI-ha': [0.1422,-0.4243,0.1426,],\
            'HDI-ha4': [0.0236,-0.2464,0.0402,],\
            'HDI-ha8': [0.0132,-0.2393,0.0267,],\
            'HDI-R': [0.0156,-0.1437,0.0024,],\
            'HDI-ha12': [0.0146,-0.2624,0.0401,],\
            'HDI-ha16': [0.0178,-0.2849,0.0547,],\
            'WFC-r': [0.0002,-0.0123,0.0027,],\
            'WFC-ha197': [0.1014,-0.3588,0.1075,],\
            'WFC-ha227': [0.0138,-0.2456,0.0303,],\
            'panstarrs-g': [0.0000,1.0000,0.0000,],\
            'panstarrs-r': [0.0000,0.0000,0.0000,],
            }
            
        # fit using stars with -0.2 < g-r < 1.2
        
        filter_trans_dict = {
            'BOK90prime-r': [0.0092,-0.1027,0.0119,],\
            'BOK90prime-ha4': [0.0389,-0.2758,0.0491,],\
            'MOS-r': [0.0027,-0.0356,0.0027,],\
            'MOS-ha4': [0.0389,-0.2758,0.0491,],\
            'MOS-ha8': [0.0178,-0.2448,0.0215,],\
            'MOS-ha12': [0.0173,-0.2683,0.0435,],\
            'MOS-R': [0.0191,-0.1550,-0.0167,],\
            'MOS-ha16': [0.0203,-0.2996,0.0739,],\
            'HDI-r': [0.0011,-0.0080,-0.0024,],\
            'HDI-ha': [0.1562,-0.4565,0.1524,],\
            'HDI-ha4': [0.0287,-0.2559,0.0392,],\
            'HDI-ha8': [0.0161,-0.2390,0.0165,],\
            'HDI-R': [0.0165,-0.1340,-0.0179,],\
            'HDI-ha12': [0.0174,-0.2736,0.0497,],\
            'HDI-ha16': [0.0215,-0.3107,0.0864,],\
            'WFC-r': [0.0002,-0.0125,0.0029,],\
            'WFC-ha197': [0.1152,-0.3910,0.1183,],\
            'WFC-ha227': [0.0161,-0.2427,0.0176,],\
            'panstarrs-g': [0.0000,1.0000,0.0000,],\
            'panstarrs-r': [0.0000,0.0000,0.0000,]
            }



        instrument_dict = {'b':'BOK90prime', 'm':'MOS','h':'HDI','i':'WFC'}

        filter_key = f"{instrument_dict[self.instrument]}-{self.filter}"

        ###################################################
        # PANSTARRS magnitudes
        ###################################################        
        PS1_r = self.pan['rmag']
        PS1_g = self.pan['gmag']
        PS1_i = self.pan['imag']        
        self.pan_gr_color = self.pan['gmag'] - self.pan['rmag']

        try:
            pcoeff = filter_trans_dict[filter_key]     
            self.R = PS1_r + pcoeff[0]*(PS1_g-PS1_r)**2 + pcoeff[1]*(PS1_g-PS1_r) + pcoeff[2]
        except KeyError:
            print(f"ruh - roh!  did not find the panstarrs color transformation for {filter_key}!!!")
            print("setting instrumental r mag to panstarrs r mag")
            print()
            self.R = self.pan['rmag']
        
    def color_correct_panstarrs_matteo(self):
        """
        NOT USED!!!

        correcting panstarrs magnitudes into the observed filter systems using conversions from M. Fossati  

        Best fit quadratic Intha - PS1_r = 0.0182*(PS1_g-PS1_r)^2 + -0.2662*(PS1_g-PS1_r) + 0.0774

        Best fit quadratic Ha4 - PS1_r = 0.0016*(PS1_g-PS1_r)^2 + -0.2134*(PS1_g-PS1_r) + 0.0168

        Best fit quadratic KPSr - PS1_r = 0.0084*(PS1_g-PS1_r)^2 + -0.0420*(PS1_g-PS1_r) + 0.0036

        Best fit quadratic KPHr - PS1_r = 0.0170*(PS1_g-PS1_r)^2 + -0.1864*(PS1_g-PS1_r) + 0.0213

        Best fit quadratic INTSr - PS1_r = 0.0023*(PS1_g-PS1_r)^2 + -0.0122*(PS1_g-PS1_r) + 0.0003


        """
        PS1_r = self.pan['rmag']
        PS1_g = self.pan['gmag']
        PS1_i = self.pan['imag']        
        self.pan_gr_color = self.pan['gmag'] - self.pan['rmag']        
        if self.filter == 'R' and ((self.instrument == 'h') | (self.instrument == 'm')): # this should be the only observations through an R filter
            print("correcting color for R filter at KPNO")            
            ###################################
            # Calculate Johnson R
            # from http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
            ###################################
            #self.R = self.pan['rmag'] + (-0.153)*(self.pan['rmag']-self.pan['imag']) - 0.117
            ###################################
            # Other transformations from 
            # https://arxiv.org/pdf/1706.06147.pdf
            # R - r = C0 + C1 x (r-i)  (-0.166, -0.275)
            # R - r = C0 + C1 x (g-r)  (-0.142, -0.166)
            ###################################
            #
            #if self.useri:
            #    self.R = self.pan['rmag'] + (-0.166)*(self.pan['rmag']-self.pan['imag']) - 0.275
            #else:
            #    self.R = self.pan['rmag'] + (-0.142)*(self.pan['gmag']-self.pan['rmag']) - 0.142

            # from Matteo Fossati
            #Best fit quadratic KPHr - PS1_r = 0.0170*(PS1_g-PS1_r)^2 + -0.1864*(PS1_g-PS1_r) + 0.0213
            self.R = PS1_r + 0.0170*(PS1_g-PS1_r)**2 + -0.1864*(PS1_g-PS1_r) + 0.0213

            # Aug 2025 - testing Kostov & Bonev linear transformation
            #self.R = self.pan['rmag'] + (-0.142)*(self.pan['gmag']-self.pan['rmag']) - 0.142
            
            # Aug 2025
            # testing ZP for HDI without using color transformation
            # trying to track down the offset we are seeing between HDI and other instruments
            # in the sense that the HDI R mags are systematically fainter

            #self.R = PS1_r            

        elif self.filter == 'r' and self.instrument == 'i':
            print("correcting color for r filter at INT")                        
            #self.R = self.pan['rmag']
            #Best fit quadratic INTSr - PS1_r = 0.0023*(PS1_g-PS1_r)^2 + -0.0122*(PS1_g-PS1_r) + 0.0003
            self.R = PS1_r + 0.0023*(PS1_g-PS1_r)**2 + -0.0122*(PS1_g-PS1_r) + 0.0003            
        # which filter is the bok telescope using?
        elif self.filter == 'r' and self.instrument == 'b':
            print("correcting color for r filter at BOK")                        
            #self.R = self.pan['rmag']
            #Best fit quadratic KPSr - PS1_r = 0.0084*(PS1_g-PS1_r)^2 + -0.0420*(PS1_g-PS1_r) + 0.0036
            self.R = PS1_r + 0.0084*(PS1_g-PS1_r)**2 + -0.0420*(PS1_g-PS1_r) + 0.0036

            # use Legacy calibration
            # https://www.legacysurvey.org/dr10/description/
            # r90Prime = rPS+0.00110−0.06875(g−i)+0.02480(g−i)^2−0.00855(g−i)^3
            ##
            ## THIS DID NOT WORK - GOING BACK TO MATTEOS TRANFORMATION
            ##
            #self.R = PS1_r  +0.00110 - 0.6875*(PS1_g-PS1_i) + 0.0248*(PS1_g-PS1_i)**2 + -0.00855*(PS1_g-PS1_i)**3  
        # this is the kpno r filter
        elif self.filter == 'r' and self.instrument == 'h':
            print("correcting color for r filter at KPNO")            
            #Best fit quadratic KPSr - PS1_r = 0.0084*(PS1_g-PS1_r)^2 + -0.0420*(PS1_g-PS1_r) + 0.0036

            
            self.R = self.pan['rmag']
            self.R = PS1_r + 0.0084*(PS1_g-PS1_r)**2 + -0.0420*(PS1_g-PS1_r) + 0.0036



        # halpha filters
        elif self.filter == 'ha' and self.instrument == 'i':
            print("correcting color for halpha filter at INT")
            #Best fit quadratic Intha - PS1_r = 0.0182*(PS1_g-PS1_r)^2 + -0.2662*(PS1_g-PS1_r) + 0.0774
            self.R = PS1_r + 0.0182*(PS1_g-PS1_r)**2 + -0.2662*(PS1_g-PS1_r) + 0.0774


            #self.R = self.pan['rmag']
        # bok is using the kpno halpha+4nm filter, so use the same correction for these
        elif self.filter == 'ha' and ((self.instrument == 'b') | (self.instrument == 'h') | (self.instrument == 'm')) :
            print("correcting color for ha filter at KPNO")                        
            #Best fit quadratic Ha4 - PS1_r = 0.0016*(PS1_g-PS1_r)^2 + -0.2134*(PS1_g-PS1_r) + 0.0168
            #self.R = self.pan['rmag']
            self.R = PS1_r + 0.0016*(PS1_g-PS1_r)**2 + -0.2134*(PS1_g-PS1_r) + 0.0168
        else:
            print("ruh - roh!  did not find the panstarrs color transformation!!!")
            print("setting instrumental r mag to panstarrs r mag")
            print()
            self.R = self.pan['rmag']
    def plot_se_pan(self):
        plt.figure()
        x = self.R[self.matchflag]
        y = self.matchedarray1['MAG_APER'][:,self.naper][self.matchflag]
        plt.plot(x,y,'bo')
        plt.xlabel('Corrected PanSTARRS Magnitude')
        plt.ylabel('SE Magnitude')

    def plot_se_pan_positions(self):
        """plot positions of SE and panstarrs catalogs  """
        plt.figure(figsize=(10,10))
        plt.plot(self.pan['RAJ2000'],self.pan['DEJ2000'],'bo',markersize=10,mfc='None',label='PanSTARRS')
        
        plt.plot(self.secat['ALPHA_J2000'],self.secat['DELTA_J2000'],'r.',label='SE')
        plt.legend()
        plt.gca().invert_xaxis()
        #print(f"number of matched sources = {np.sum(zp.matchflag)}")

        # add circle

        circle1 = plt.Circle((self.centerRA, self.centerDEC), self.radius, color='c',alpha=.2)
        plt.gca().add_patch(circle1)

        plt.savefig('plots/'+self.plotprefix.replace('.fits','')+'se-pan-positions.png')        
        
    def plot_fitresults(self, x, y, yerr=None, polyfit_results = [0,0],color=None,ymin_fixed=-0.12,ymax_fixed=0.12):
        # plot best-fit results
        #print(f"inside plot_fitresults, len(x) = {len(x)}, len(color)={len(color)}")
        yfit = np.polyval(polyfit_results,x)
        residual = (yfit - y)
        plt.figure(figsize=(8,8))

        self.fitmad = MAD2(residual)
        self.fitstd = np.std(residual)
        s = ' (MAD =%.4f, std=%.4f, N=%d)'%(MAD2(residual), np.std(residual),len(residual))
        
        
        plt.subplot(2,1,1)
        if len(yerr) < len(y):
            
            if color is not None:
                plt.scatter(x,residual,s=30,label=s,c=color)
                plt.colorbar(label='g-r')                
            else:
                plt.plot(x,residual, 'ko',label=s)
            ymin, ymax = plt.ylim()
                
        else:
            #print(f"inside plot_fitresults, len(x) = {len(x)}, len(color)={len(color)}")
            # plot datapoints and get ymin and ymax, then use this to get limits for plot with errorbars
            if color is not None:
                plt.scatter(x,y,s=30,label=s,c=color,zorder=10)
                plt.colorbar(label='g-r')
            else:
                plt.scatter(x,y,s=30,label=s,zorder=10)
            ymin,ymax = plt.ylim()
            plt.errorbar(x,y,yerr=yerr,fmt='None',ecolor='0.5',label='SE MAG',alpha=.4,zorder=1)
        plt.ylim(ymin,ymax) # let ymin, ymax
        plt.xlabel('Pan-STARRS Corrected',fontsize=16)
        plt.ylabel('SE MAG',fontsize=16)
        xl = np.linspace(14,19,10)
        yl = np.polyval(polyfit_results,xl)
        s = 'fit: y = %.2f PAN + %.2f'%(polyfit_results[0],polyfit_results[1])
        plt.plot(xl,yl,'k--',label=s)
        plt.legend()
        plt.title(self.plotprefix+s)
        
        plt.subplot(2,1,2)
        s = 'MAD = %.4f'%(MAD2(residual))
        if len(yerr) < len(y):

            # color by panstarrs g-r color to check for residual color terms
            if color is not None:
                plt.scatter(x,residual,s=30,label=s,c=color)
            else:
                plt.plot(x,residual, 'ko',label=s)
            ymin,ymax = plt.ylim()
        else:
            if color is not None:
                plt.scatter(x,residual,s=30,label=s,c=color,zorder=10)
                plt.colorbar(label='g-r')
            else:
                plt.scatter(x,residual,s=30,label=s,zorder=10)
            ymin,ymax = plt.ylim()
            plt.errorbar(x,residual,yerr=yerr,fmt='None',ecolor='0.5',label='SE MAG '+s,alpha=.4,zorder=1)
        plt.ylim(ymin_fixed,ymax_fixed)
        plt.xlabel('Pan-STARRS Corrected',fontsize=16)
        plt.ylabel('YFIT - SE R-band MAG',fontsize=16)
        plt.legend()
        plt.axhline(y=0,color='r')
        plt.savefig('plots/'+self.plotprefix.replace('.fits','')+'se-pan-flux.png')
        #plt.close()
    def fitzp(self,plotall=False):
        ###################################
        # Solve for the zeropoint
        ###################################

        # plot Pan-STARRS r mag on x axis, observed R-mag on y axis

        if plotall:
            flag = self.fitflag
            plt.figure(figsize=(6,4))
            plt.title(self.plotprefix)
            plt.plot(self.R[flag],self.R[flag]-self.matchedarray1['MAG_AUTO'][flag],'bo')
            plt.errorbar(self.R[flag],self.R[flag]-self.matchedarray1['MAG_AUTO'][flag],xerr= self.pan['e_rmag'][flag],yerr=self.matchedarray1['MAGERR_AUTO'][flag],fmt='none')
            plt.plot(self.R[flag],self.R[flag]-self.matchedarray1['MAG_BEST'][flag],'ro',label='MAG_BEST')
            plt.plot(self.R[flag],self.R[flag]-self.matchedarray1['MAG_PETRO'][flag],'go',label='MAG_PETRO')
            plt.plot(self.R[flag],self.R[flag]-self.matchedarray1['MAG_APER'][:,self.naper][flag],'ko',label='MAG_APER')
            plt.xlabel('Pan-STARRS r',fontsize=16)
            plt.ylabel('PANSTARRS - SE_MAG',fontsize=16)
            plt.legend()
            xl = np.linspace(14,17,10)
            yl = 0*xl
            #yl = np.polyval(c,xl)
            #plt.plot(xl,yl,'k--')
            plt.savefig('getzp-fig2.png')
            #plt.plot(xl,1.2*yl,'k:')

        flag = self.fitflag
        residual = np.zeros(len(flag))
        ####################################
        ## had been dividing by yfit, but that doesn't make sense
        ## want residual to be in magnitudes
        ## removing yfit normalization
        ####################################

        # fixed radii apertures: [:,0] = 3 pix, [:,1] = 5 pix, [:,2] = 7 pixels

        # TODONE - incorporate panstarrs errors in the yerr
        yerrpan = self.pan['e_rmag'][flag]
        if self.mag == 0: # this is the default magnitude
            if self.verbose:
                print('Using Aperture Magnitudes')
            y = self.matchedarray1['MAG_APER'][:,self.naper][flag]
            yerr = np.sqrt(self.matchedarray1['MAGERR_APER'][:,self.naper][flag]**2 + yerrpan**2)
        elif self.mag == 1:
            if self.verbose:
                print('Using MAG_BEST')
            y = self.matchedarray1['MAG_BEST'][flag]
            yerr = np.sqrt(self.matchedarray1['MAGERR_BEST'][flag]**2 + yerrpan**2)
        elif self.mag == 2:
            if self.verbose:
                print('Using MAG_PETRO')
            y = self.matchedarray1['MAG_PETRO'][flag]
            
            yerr = np.sqrt(self.matchedarray1['MAGERR_PETRO'][flag]**2 + yerrpan**2)


        # start fitting procedure
        flag = self.fitflag
        self.bestc = np.array([0,0],'f')
        delta = 100.     
        x = self.R[flag] # expected mag from panstarrs
        color = self.pan_gr_color[flag]
        #print(f"len(x) = {len(x)}, len(color)= {len(color)}")

        useAstropyModel = args.useastropy

        if useAstropyModel:
            print("using astropy model...")

            # get initial fit from numpy inside fitline_sigma_clipping function
            fit2, mask2,uncertainty = fitline_sigma_clipping(x,y, yerr=yerr,niter=5)#,sigma = yerr)
            slope = 1

            yplot = self.matchedarray1['MAG_APER'][:,self.naper][self.fitflag]
            magfit = fit2(self.R[self.fitflag])
            residual_all = 10.**((magfit - yplot)/2.5)        
            self.residual_all = residual_all

            residual_match = residual_all[~mask2]
            x = x[~mask2]
            y = y[~mask2]
            yerr = yerr[~mask2]
            color = color[~mask2]

            #print("checking length of x,y, residual_all: ", len(x),len(y), len(self.residual_all)) 
            self.zpcovar = [[0,0],[0,uncertainty[0]**2]]
            self.zperr = uncertainty[0]
            self.zp = fit2.intercept.value
            #self.fitflag[self.fitflag] = self.fitflag[self.fitflag] & mask2
            self.bestc = [1,self.zp]
            #self.plot_fitresults(x,y,yerr=yerr,polyfit_results = self.bestc,color=residual_match)
            self.plot_fitresults(x,y,yerr=yerr,polyfit_results = self.bestc,color=color)
            
        else:
            
            while delta > 1.e-3:
                #c = np.polyfit(x,y,1)
                t = curve_fit(zpfunc,x,y,sigma = yerr)
                # convert to format expected from polyfit
                c = np.array([1.,t[0][0]])
                if self.verbose:
                    print('number of points retained = ',np.sum(flag))
                yfit = np.polyval(c,x)
                residual = (yfit - y)

                if plotall:
                    self.plot_fitresults(x,y,yerr=yerr,polyfit_results = c,color=color)


                # check for convergence
                if self.verbose:
                    print('new ZP = {:.3f}, previous ZP = {:.3f}'.format(self.bestc[1],c[1]))
                delta = abs(self.bestc[1] - c[1])
                self.bestc = c
                MAD = mad_std(residual)#1.48*np.median(abs(residual - np.median(residual)))
                #clip_flag = sigma_clip(self.zim,sigma=3,maxiters=10,masked=True)            
                flag =  (abs(residual - np.median(residual)) < self.nsigma*MAD)
                if sum(flag) < 2:
                    print(f'WARNING: ONLY ONE DATA POINT LEFT FOR {self.image}')
                    self.x = x
                    self.y = y
                    self.residual = residual
                    sys.exit()
                #flag =  (abs(residual) < self.nsigma*np.std(residual))
                self.x = x
                self.y = y
                self.residual =residual
                x = x[flag]
                y = y[flag]
                yerr = yerr[flag]
                color = color[flag]
            

            ###################################
            ##  show histogram of residuals
            ###################################

            yplot = self.matchedarray1['MAG_APER'][:,self.naper][self.fitflag]
            magfit = np.polyval(self.bestc,self.R[self.fitflag])
            residual_all = 10.**((magfit - yplot)/2.5)        
            self.residual_all = residual_all

            self.zpcovar = t[1]
            self.zperr = np.sqrt(self.zpcovar[0][0])
            self.zp = self.bestc[1]
            self.plot_fitresults(x,y,yerr=yerr,polyfit_results = self.bestc,color=color)            
        s = 'residual (mean,std) = %.3f +/- %.3f'%(np.mean(residual_all),np.std(residual_all))
        if self.verbose:
            print(s)
        if plotall:
            plt.figure()            
            crap = plt.hist(residual_all,bins=np.linspace(.8,1.8,20))
            plt.text(0.05,.85,s,horizontalalignment='left',transform=plt.gca().transAxes)
            plt.xlabel('Residuals')
            plt.savefig('plots/'+self.plotprefix.replace(".fits","")+'getzp-residual-hist.png')



        ###################################
        # Show location of residuals
        ###################################
        '''
        # this plots locations of all sources, not just the ones that 
        # are used in the ZPfitting
        # 
        plt.figure(figsize=(6,4))
        plt.title(self.image)
        yplot2 = self.matchedarray1['MAG_APER'][:,self.naper]
        magfit2 = np.polyval(self.bestc,self.R)
        residual_all2 = 10.**((magfit2 - yplot2)/2.5)

        plt.scatter(self.matchedarray1['X_IMAGE'],self.matchedarray1['Y_IMAGE'],c = (residual_all2),vmin=v1,vmax=v2,s=15)
        cb=plt.colorbar()
        cb.set_label('f-WFC/f-pan')        
        plt.savefig('getzp-position-residuals-all-fig1.png')
        '''
        plt.figure(figsize=(6,4))
        
        s = ' (mean,std,MAD = {:.3f},{:.3f},{:.3f})'.format(np.mean(residual_all),np.std(residual_all),MAD2(residual_all))
        #s = str(MAD(residual_all))
        plt.title(self.plotprefix+s)
        self.residual_allx = self.matchedarray1['X_IMAGE'][self.fitflag]
        self.residual_ally = self.matchedarray1['Y_IMAGE'][self.fitflag]

        # debugging some wacky BOK results - want to see full range of residuals
        # so letting vmin/vmax get set automatically
        # issue is with astropy sigma clipping - ugh!!!
        # returning to normal...
        plt.scatter(self.matchedarray1['X_IMAGE'][self.fitflag],self.matchedarray1['Y_IMAGE'][self.fitflag],c = (residual_all),vmin=v1,vmax=v2,s=15)
        #plt.scatter(self.matchedarray1['X_IMAGE'][self.fitflag],self.matchedarray1['Y_IMAGE'][self.fitflag],c = (residual_all),vmin=.5,vmax=2.5,s=15)
        cb=plt.colorbar()
        cb.set_label('f-meas/f-pan')
        plt.savefig('plots/'+self.plotprefix.replace(".fits","")+'getzp-xyresidual-fitted.png')

        self.x = x
        self.y = y
        self.yerr = yerr

        

    def check_90prime_ccds(self):
        image_med = np.ma.median(self.residual_all)
        hdu = fits.open(self.image)                
        NAXIS1 = self.header['NAXIS1']
        NAXIS2 = self.header['NAXIS2']
        quad = 1

        for ix in range(2):
            xmin = 0 + NAXIS1//2*ix
            xmax = NAXIS1//2 * (ix +1)
            
            for iy in range(2):
                ymin = 0 + NAXIS2//2*iy
                ymax = NAXIS2//2 * (iy+1)
                #print(xmin,xmax,ymin,ymax)
                flag = (self.residual_allx > xmin) & (self.residual_allx < xmax) &\
                    (self.residual_ally > ymin) & (self.residual_ally < ymax) 
                amp_med = np.ma.median(self.residual_all[flag])
                
                amp_scale = image_med/amp_med
                
                print(f"\tmedian and scale for quadrant {quad} = {amp_med:.3f} {amp_scale:.3f}")
                # scale the data
                # scale the image and ivar accordingly            
                hdu[0].data[ymin:ymax,xmin:xmax] = amp_scale*hdu[0].data[ymin:ymax,xmin:xmax]

                # what is the correct way to scale the weights?  same as image or inverse???
                # in the weight image, high numbers are good
                quad += 1
        hdu.writeto(self.image, overwrite=True)

    def check_int_ccds(self, trim=100, min_stars=5, max_scale_dev=None):
        """
        Renormalize INT/WFC CCD regions using per-star residuals.

        INT chip layout in mosaic coordinates:
          CCD3: x=   1:4300, y=4560:6420
          CCD4: x=   1:4300, y=2430:4280
          CCD1: x=   1:4300, y=   1:2300
          CCD2: x=4420:6450, y=   1:4420   (sideways chip)

        Parameters
        ----------
        trim : int, optional
            Pixels to trim from each chip boundary when selecting stars to
            estimate the chip normalization. The full chip region is still
            scaled after the factor is derived.
        min_stars : int, optional
            Minimum number of stars required on a chip to derive a correction.
        max_scale_dev : float or None, optional
            Optional safety limit. If set, skip applying a correction when
            abs(scale - 1) > max_scale_dev.
            Example: max_scale_dev=0.15 skips scales outside [0.85, 1.15].
        """

        # Full chip footprints (used when applying the correction)
        chip_boxes = [
            ("CCD3", (1,    4300, 4560, 6420)),
            ("CCD4", (1,    4300, 2430, 4280)),
            ("CCD1", (1,    4300, 1,    2300)),
            ("CCD2", (4420, 6450, 1,    4420)),
        ]

        # Trimmed chip footprints (used when measuring the correction)
        trimmed_boxes = []
        for name, (xmin, xmax, ymin, ymax) in chip_boxes:
            x0 = xmin + trim
            x1 = xmax - trim
            y0 = ymin + trim
            y1 = ymax - trim

            if (x1 <= x0) or (y1 <= y0):
                print(f"{name}: trim={trim} too large; skipping chip")
                continue

            trimmed_boxes.append((name, (x0, x1, y0, y1)))

        if len(trimmed_boxes) == 0:
            print("No valid INT chip regions after trimming; skipping correction")
            return

        # Build mask of stars that land on any trimmed chip region
        on_chip = np.zeros(len(self.residual_all), dtype=bool)
        for name, (x0, x1, y0, y1) in trimmed_boxes:
            on_chip |= (
                (self.residual_allx >= x0) & (self.residual_allx < x1) &
                (self.residual_ally >= y0) & (self.residual_ally < y1)
            )

        n_on_chip = np.sum(on_chip)
        if n_on_chip < min_stars:
            print(f"Only {n_on_chip} stars on trimmed INT chip regions; skipping correction")
            return

        image_med = clipped_median(self.residual_all[on_chip])
        if not np.isfinite(image_med) or image_med == 0:
            print(f"Invalid global INT residual median ({image_med}); skipping correction")
            return

        print(f"Global median residual over INT chips = {image_med:.4f}")

        hdu = fits.open(self.image)
        applied_scales = {}

        for name, (xmin, xmax, ymin, ymax) in chip_boxes:
            # matching trimmed box for measurement
            tmatch = [box for n, box in trimmed_boxes if n == name]
            if len(tmatch) == 0:
                print(f"\t{name}: no valid trimmed region; skipping")
                continue
            x0, x1, y0, y1 = tmatch[0]

            flag = (
                (self.residual_allx >= x0) & (self.residual_allx < x1) &
                (self.residual_ally >= y0) & (self.residual_ally < y1)
            )

            nstars = np.sum(flag)
            if nstars < min_stars:
                print(f"\t{name}: only {nstars} stars; skipping")
                continue

            chip_med = clipped_median(self.residual_all[flag])
            if not np.isfinite(chip_med) or chip_med == 0:
                print(f"\t{name}: invalid median {chip_med}; skipping")
                continue

            chip_scale = image_med / chip_med

            if (max_scale_dev is not None) and (abs(chip_scale - 1.0) > max_scale_dev):
                print(f"\t{name}: scale={chip_scale:.4f} exceeds limit; skipping")
                continue

            print(f"\t{name}: nstars={nstars:3d} median={chip_med:.4f} scale={chip_scale:.4f}")

            # Apply correction to the full nominal chip footprint
            hdu[0].data[ymin:ymax, xmin:xmax] *= chip_scale
            applied_scales[name] = chip_scale

        # Optional bookkeeping in the header
        hdr = hdu[0].header
        hdr["INTFIX"] = (True, "INT CCD normalization applied")
        for name in ["CCD1", "CCD2", "CCD3", "CCD4"]:
            key = f"I{name[-1]}SCL"
            if name in applied_scales:
                hdr[key] = (float(applied_scales[name]), f"{name} scale factor")
            else:
                hdr[key] = (-999.0, f"{name} scale factor not applied")

        hdu.writeto(self.image, overwrite=True)
        hdu.close()


    def update_header(self):
        #print('working on this')
        # add best-fit ZP to image header
        im, header = fits.getdata(self.image,header=True)
        zperr = self.zperr         
        # or convert vega zp to AB
        if self.filter == 'R':
            # conversion from Blanton+2007
            # http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
            # header.set('PHOTZP',float('{:.3f}'.format(-1.*self.bestc[1]+.21)))
            # changed this to write out phot zp in AB system for ALL filters
            header.set('PHOTZP',float('{:.3f}'.format(-1.*self.zp)),'ZP fit from panstarrs AB mag')

            header.set('LAMB_um',float(.6442))

        else:
            header.set('PHOTZP',float('{:.3f}'.format(-1.*self.zp)),'ZP fit from panstarrs AB mag')
        print(f"PHOTZPER = {zperr}")
        try:
            header.set('PHOTZPER',float(f'{zperr:.4f}'),'ZP fit err from covariance')
        except ValueError:
            header.set('PHOTZPER',f'{zperr}','ZP fit err from covariance')            
        try:
            header.set('PHOTZMAD',float(f'{self.fitmad:.3f}'),'ZP fit MAD of residuals in mag')
        except ValueError:
            print("problem adding to header: MAD of the fit")            
        try:
            header.set('PHOTZSTD',float(f'{self.fitstd:.3f}'),'ZP fit STD of residuals in mag')
        except ValueError:
            print("problem adding to header: STD of the fit")            
            
        header.set('PHOTSYS','AB')
        header.set('FLUXZPJY',float(3631))
        # add flatten number to header
        npreviousflat = 0
        suffixes = ['','1','2']
        for s in suffixes:
            try:
                temp = header['ZPNFLAT'+s]
                npreviousflat += 1
            except KeyError:
                pass
        # add a new header keyword if the image has already been flattened
        # this could happen with renaming coadds to VFS format and then rerunning getzp
        #header.set('ZPNFLAT'+suffixes[npreviousflat],int(self.flatten),'getzp --flatten number')
        header.set('ZPNFLAT',int(self.flatten),'getzp --flatten number')        
        fits.writeto(self.image, im, header, overwrite=True)


    def fit_residual_surface_old(self, norder=2, suffix=None):
        mag_meas = self.matchedarray1['MAG_APER'][:, self.naper][self.fitflag]
        mag_pan = np.polyval(self.bestc, self.R[self.fitflag])
        flux_ratio = 10.**((mag_pan - mag_meas) / 2.5)

        self.xim = self.matchedarray1['X_IMAGE'][self.fitflag]
        self.yim = self.matchedarray1['Y_IMAGE'][self.fitflag]
        self.zim = flux_ratio

        self.imagedata = fits.getdata(self.image)
        ny, nx = self.imagedata.shape

        clip_flag = sigma_clip(self.zim, sigma=3, maxiters=10, masked=True)
        good = ~clip_flag.mask

        if self.spline:
            # depends on fitspline2d convention, but ideally this should also
            # evaluate on the true detector grid
            self.zz = fitspline2d(
                self.xim[good], self.yim[good], self.zim[good],
                nx, ny, order=self.spline_order, s=self.spline_smooth
            )
        else:
            m = polyfit2d(self.xim[good], self.yim[good], self.zim[good], order=norder)

            # evaluate on true detector coordinates
            yy, xx = np.indices((ny, nx))
            self.zz = polyval2d(xx, yy, m)

        print("zim stats:", np.nanmin(self.zim[good]), np.nanmedian(self.zim[good]), np.nanmax(self.zim[good]))
        print("zz stats :", np.nanmin(self.zz), np.nanmedian(self.zz), np.nanmax(self.zz))
        good = ~clip_flag.mask
        print("N stars used:", np.sum(good))
        print("zim median/std/MAD:",
                  np.nanmedian(self.zim[good]),
                  np.nanstd(self.zim[good]),
                  MAD2(self.zim[good]))
        print("zz min/median/max:",
                  np.nanmin(self.zz),
                  np.nanmedian(self.zz),
                  np.nanmax(self.zz))
    def fit_residual_surface(self, norder=2, suffix=None):

        self.xim = self.residual_allx
        self.yim = self.residual_ally
        self.zim = self.residual_all

        self.imagedata = fits.getdata(self.image)
        ny, nx = self.imagedata.shape

        clip_flag = sigma_clip(self.zim, sigma=3, maxiters=10, masked=True)
        good = ~clip_flag.mask

        if self.spline:
            self.zz = fitspline2d(
                self.xim[good], self.yim[good], self.zim[good],
                nx, ny, order=self.spline_order, s=self.spline_smooth
            )
        else:
            m = polyfit2d(self.xim[good], self.yim[good], self.zim[good], order=norder)
            yy, xx = np.indices((ny, nx))
            self.zz = polyval2d(xx, yy, m)

        print("N stars used:", np.sum(good))
        print("zim median/std/MAD:",
              np.nanmedian(self.zim[good]),
              np.nanstd(self.zim[good]),
              MAD2(self.zim[good]))
        print("zz min/median/max:",
              np.nanmin(self.zz),
              np.nanmedian(self.zz),
              np.nanmax(self.zz))
    
        print("plotting results of background fitting...")
        plt.figure(dpi=150)
        n = 10 # downsampling to speed up 
        plt.imshow(self.zz[::n,::n],vmin=v1,vmax=v2,origin="lower")
        #plt.imshow(self.zz,vmin=v1,vmax=v2,origin="lower")
        plt.scatter(self.xim[~clip_flag.mask]/n, self.yim[~clip_flag.mask]/n, c=self.zim[~clip_flag.mask],vmin=v1,vmax=v2,s=15)
        cb=plt.colorbar()
        cb.set_label('f-meas/f-pan')
        s = ' std (MAD) = %.4f (%.4f)'%(np.std(self.zim[~clip_flag.mask]),MAD2(self.zim[~clip_flag.mask]))
        plt.title(self.plotprefix+': n poly = '+str(norder)+s)
        #plt.show()
        plt.xlabel(f"pixel/{n}")
        plt.ylabel(f"pixel/{n}")        
        if suffix is None:
            plotname='imsurfit-'+str(norder)
        else:
            plotname='imsurfit-'+str(norder)+'-'+suffix
        plt.savefig('plots/'+self.plotprefix.replace(".fits","")+plotname+'.png')
        #plt.savefig('plots/'+self.plotprefix.replace(".fits","")+plotname+'.pdf')
        print("... done plotting results.\n")

    def fit_residual_surface_old_old(self,norder=2,suffix=None):
        """
        for INT data, first pass:
        - create an image of the fit the residuals WRT panstarrs,
        - use SE to fit the background to the image (or some other polynomial fit?)
        - normalize the background image
        - divide this into the science frame
        - resolve for the photometric zp

        """
        # get difference between measured and predicted fluxes
        mag_meas = self.matchedarray1['MAG_APER'][:,self.naper][self.fitflag]
        mag_pan = np.polyval(self.bestc,self.R[self.fitflag])
        flux_ratio = 10.**((mag_pan - mag_meas)/2.5)
        
        # create an image
        self.xim = self.matchedarray1['X_IMAGE'][self.fitflag]
        self.yim = self.matchedarray1['Y_IMAGE'][self.fitflag]        

        self.zim = flux_ratio
        self.imagedata = fits.getdata(self.image)
        # fill in where there is no coverage        
        weight_name = self.image.split('.fits')[0]+'.weight.fits'
        self.weightdata = fits.getdata(weight_name)
        # this is not used
        self.nodata =  self.weightdata == 0


        if self.instrument == 'i':

            ##########################################################
            # This is some fine tuning for INT data, but not actually using this?
            #
            # center of geometric distortion ~= (3500, 2400)
            # top chip has y > 4300
            # blank area starts at x > 4300
            # end of image at x = 6300, y=6400
            # so drop some of top left points into top right
            # 
            flip_data = (self.xim < (6300-3500)) & (self.yim > 4300)

            # add points from the top left corner into blank corner
            # looks like I'm not actually using this though
            nrandom = 40
            fakex = 2*3500 - self.xim[flip_data]
            fakey = self.yim[flip_data]
            fakez = self.zim[flip_data]                

            fakey = fakey[fakex < 6300]
            fakez = fakez[fakex < 6300]                
            fakex = fakex[fakex < 6300]        

            # combine fake data with original data
            #self.xim = np.array(self.xim.tolist()+fakex.tolist())
            #self.yim = np.array(self.yim.tolist()+fakey.tolist())
            #self.zim = np.array(self.zim.tolist()+fakez.tolist())
            ##########################################################

            
        # clip data
        clip_flag = sigma_clip(self.zim,sigma=3,maxiters=10,masked=True)
        ny, nx = self.imagedata.shape
        if self.spline:
            self.zz = fitspline2d(self.xim[~clip_flag.mask],self.yim[~clip_flag.mask],self.zim[~clip_flag.mask],nx,ny,order=self.spline_order,s=self.spline_smooth)

        else:
            # Fit a 2nd order, 2d polynomial
            m = polyfit2d(self.xim[~clip_flag.mask],self.yim[~clip_flag.mask],self.zim[~clip_flag.mask],order=norder)

            # Evaluate it on a grid...        

            #xx, yy = np.meshgrid(np.linspace(self.xim.min(), self.xim.max(), nx), 
            #             np.linspace(self.yim.min(), self.yim.max(), ny))
            yy, xx = np.indices((ny, nx))
            self.zz = polyval2d(xx, yy, m)
            

            # Plot
        print("zz stats:", np.nanmin(self.zz), np.nanmedian(self.zz), np.nanmax(self.zz))
        print("zim stats:", np.nanmin(self.zim), np.nanmedian(self.zim), np.nanmax(self.zim))
        
        print("plotting results of background fitting...")
        plt.figure(dpi=150)
        n = 10 # downsampling to speed up 
        plt.imshow(self.zz[::n,::n],vmin=v1,vmax=v2,origin="lower")
        #plt.imshow(self.zz,vmin=v1,vmax=v2,origin="lower")
        plt.scatter(self.xim[~clip_flag.mask]/n, self.yim[~clip_flag.mask]/n, c=self.zim[~clip_flag.mask],vmin=v1,vmax=v2,s=15)
        cb=plt.colorbar()
        cb.set_label('f-meas/f-pan')
        s = ' std (MAD) = %.4f (%.4f)'%(np.std(self.zim[~clip_flag.mask]),MAD2(self.zim[~clip_flag.mask]))
        plt.title(self.plotprefix+': n poly = '+str(norder)+s)
        #plt.show()
        plt.xlabel(f"pixel/{n}")
        plt.ylabel(f"pixel/{n}")        
        if suffix is None:
            plotname='imsurfit-'+str(norder)
        else:
            plotname='imsurfit-'+str(norder)+'-'+suffix
        plt.savefig('plots/'+self.plotprefix.replace(".fits","")+plotname+'.png')
        #plt.savefig('plots/'+self.plotprefix.replace(".fits","")+plotname+'.pdf')
        print("... done plotting results.\n")


    def renorm_wfc(self):
        self.zz_norm = self.zz / np.nanmedian(self.zz)

        bad = ~np.isfinite(self.zz_norm) | (self.zz_norm <= 0)
        if np.any(bad):
            print(f"WARNING: replacing {np.sum(bad)} bad correction pixels with 1.0")
            self.zz_norm = np.where(bad, 1.0, self.zz_norm)

        self.imagedata, header = fits.getdata(self.image, header=True)
        self.imagedata_norm = self.imagedata / self.zz_norm
        print("DEBUG renorm_wfc input image:", self.image)
        print("DEBUG correction min/med/max:",
                  np.nanmin(self.zz_norm), np.nanmedian(self.zz_norm), np.nanmax(self.zz_norm))

        ratio = self.imagedata_norm / self.imagedata
        print("DEBUG flattened/original min/med/max:",
                  np.nanmin(ratio), np.nanmedian(ratio), np.nanmax(ratio))

        self.renorm_image = 'f' + self.image
        print("DEBUG renorm_wfc output image:", self.renorm_image)
    
        #self.renorm_image = 'f' + self.image
        fits.writeto(self.renorm_image, self.imagedata_norm, header=header, overwrite=True)

    
    def rerun_zp_fit(self):
        # change image name to flattened image
        self.image = self.renorm_image
        self.plotprefix = 'f'+self.plotprefix
        # rerun getzp, but don't download panstarrs again
        self.runse(useprevious=False)
        if self.verbose:
            print('STATUS: matching se cat to panstarrs')
            
        self.match_coords()
        #self.color_correct_panstarrs()        

        if self.verbose:
            print('STATUS: fitting zeropoint')
            
        self.fitzp()
        if self.verbose:
            print('STATUS: udating header')        
        self.update_header()
        # check to make sure the systematic residuals have been removed
        self.fit_residual_surface(suffix='round2',norder=self.norder)

    def write_pan_se_table(self):

        # join the panstarrs and SE matched tables
        outtab = hstack([self.pan,Table(self.matchedarray1)])

        # trying to fix error that we were getting saying that keyword "description" is too long
        # (>8 characters or has illegal characters)
        try:
            original_description = outtab.meta['description']
            del outtab.meta['description']
            outtab.meta['DSCRPTN'] = original_description[:70] # Truncating the value to 70 characters 
        except KeyError:
            print("WARNING: could not get description in meta data.  Maybe this is ok?")

        
        ##################################################
        # convert the SE mag columns using best fit zp
        ##################################################        
        #print("columns in pan se merged table: ", outtab.colnames)
        # magnitude colums
        magcols = ['MAG_'+i for i in ['ISO','ISOCOR','AUTO','BEST','PETRO']]
        for m in magcols:
            outtab[m] = outtab[m] +  -self.zp

        # aperture magnitude columns
        n_aper_mag = len(outtab['MAG_APER'][0])
        for i in range(n_aper_mag):
            outtab['MAG_APER'][:,i] += -self.zp 
            

        # write out combined table
        subdir = 'matched_panstarrs_se_tables'
        if not os.path.exists(subdir):
            os.mkdir(subdir)

        #print("HEY!!! self.image = ",self.image)
        #print("HEY!!! get_filebasename(self.image) = ",get_filebasename(self.image))
        fname = f"{get_filebasename(self.image)}_{self.filter}_pan_SE_tab.fits"
        #print("table name = ",fname)
        outname = os.path.join(subdir, fname)
        #outname = get_filebasename(self.image)+'_pan_SE_tab.fits'
        #print("Writing merged panstarrs - SE table as ",outname)
        # only keep stars that are used in fitting for the ZP
        outtab = outtab[self.fitflag]
        # add the color-transformed R magnitude
        c = Column(self.R[self.fitflag], name='pan2instmag')
        outtab.add_column(c)
        
        outtab.write(outname, format='fits', overwrite=True)

def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Run sextractor, get Pan-STARRS catalog, and then compute photometric ZP\n"
            "\n"
            "from within ipython:\n"
            "%run ~/github/Virgo/programs/getzp.py --image pointing031-r.coadd.fits --instrument i\n"
            "\n"
            "The y intercept is -1*ZP.\n"
            "\n"
            "x and y data can be accessed at zp.x and zp.y in case you want to make additional plots."
        ),
        formatter_class=argparse.RawTextHelpFormatter
    )

    parser.add_argument('--image', dest='image', default='test.coadd.fits',
                        help='Image for ZP calibration')
    parser.add_argument('--instrument', dest='instrument', default=None,
                        help='HDI = h, KPNO mosaic = m, INT = i, BOK 90Prime = b')
    parser.add_argument('--catalog', dest='catalog', default=None,
                        help='photometric catalog to use for bootstrapping photometry')
    parser.add_argument('--fwhm', dest='fwhm', default=None,
                        help='image FWHM in arcseconds. Default is none, then SE assumes 1.5 arcsec')
    parser.add_argument('--filter', dest='filter', default='R',
                        help='filter, options are: ha4, ha8, ha12, ha16, ha197 (INT Halpha), ha227 (INT Ha6657), r, R')
    parser.add_argument('--useri', dest='useri', default=False, action='store_true',
                        help='Use r->R transformation as a function of r-i rather than the g-r relation. g-r is the default.')
    parser.add_argument('--normbyexptime', dest='normbyexptime', default=False, action='store_true',
                        help='set this flag if the image is in ADU rather than ADU/s, and the program will then normalize by the exposure time. Note: swarp produces images in ADU/s, so this is usually not necessary if using coadds from swarp.')
    parser.add_argument('--mag', dest='mag', default='0', choices=['0', '1', '2'],
                        help='select SE magnitude to use when solving for ZP. 0=MAG_APER, 1=MAG_BEST, 2=MAG_PETRO. Default is MAG_APER')
    parser.add_argument('--minmag', dest='minmag', default=15,
                        help='bright magnitude cutoff; default is 15.')
    parser.add_argument('--naper', dest='naper', default=5,
                        help='select fixed aperture magnitude. 0=10pix, 1=12pix, 2=15pix, 3=20pix, 4=25pix, 5=30pix. Default is 5 (30 pixel diameter)')
    parser.add_argument('--nsigma', dest='nsigma', default=3.5,
                        help='number of std to use in iterative rejection of ZP fitting. default is 3.5')
    parser.add_argument('--d', dest='d', default='~/github/HalphaImaging/astromatic/',
                        help='Locates path of default config files. Default is ~/github/HalphaImaging/astromatic')
    parser.add_argument('--fit', dest='fitonly', default=False, action='store_true',
                        help='Do not run SE or download catalog. just redo fitting.')
    parser.add_argument('--flatten', dest='flatten', default='0', choices=['0', '1', '2'],
                        help='Number of times to run flattening process to try to remove vignetting/illumination patterns. The default is zero. Options are [0,1,2]. This is needed for INT data from 2019. HDI does not show this effect, and INT data from 2022 does not seem to show it either.')
    parser.add_argument('--useastropy', dest='useastropy', default=False, action='store_true',
                        help='Use astropy to fit the ZP line with sigma clipping. Default is False, which then uses my own iterative fitting.')
    parser.add_argument('--spline', dest='spline', default=False, action='store_true',
                        help='Fit surface with a spline rather than a 2d polynomial')
    parser.add_argument('--spline_order', dest='spline_order', default=3,
                        help='order of spline. default is 3.')
    parser.add_argument('--spline_smooth', dest='spline_smooth', default=1000,
                        help='smoothing for spline. default is 1000.')
    parser.add_argument('--norder', dest='norder', default='2', choices=['0', '1', '2', '3', '4'],
                        help='degree of polynomial to fit to overall background. default is 2.')
    parser.add_argument('--verbose', dest='verbose', default=False, action='store_true',
                        help='print extra debug/status statements')
    parser.add_argument('--getrefcatonly', dest='getrefcatonly', default=False, action='store_true',
                        help='download reference PANSTARRS catalog only. use this before running with slurm')
    parser.add_argument('--fixint', dest='fixint', default=False, action='store_true',
                        help='fix offset b/w INT ccds. default is False.')
    parser.add_argument('--fixbok', dest='nofixbok', default=False, action='store_true',
                        help='do NOT fix offset b/w bok amps')

    return parser


def main(raw_args=None, write_table=False):
    args = build_parser().parse_args(raw_args)


    zp = getzp(args)
    zp.getzp()

    print('ZP = {:.3f} +/- {:.3f}, {}'.format(-1 * zp.zp, zp.zperr, zp.image))

    if write_table:
        print("writing out the merged panstarrs + SE table")
        zp.write_pan_se_table()

    return zp, -1 * zp.zp, zp.zperr


if __name__ == '__main__':
    start_time = time.time()
    zp, zpval, zperr = main(write_table=True)
    print(f"--- {(time.time() - start_time):.1f} seconds ---")
        

 
