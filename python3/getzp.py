#!/usr/bin/env python

'''
USAGE:

from within ipython:
%run ~/github/HalphaImaging/getzp.py --image pointing031-r.coadd.fits --instrument i --filter r

The y intercept is -1*ZP

To print the value in ipython, type:

-1*zp.bestc[1]


UPDATES:
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


from matplotlib import pyplot as plt
import matplotlib.patches as patches
import matplotlib
#matplotlib.use("Tkagg")
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.stats import sigma_clip
from astropy.table import Table

from scipy.optimize import curve_fit
try:
    from scipy.stats import median_abs_deviation as MAD2
except:
    from scipy.stats import median_absolute_deviation as MAD2
from astroquery.vizier import Vizier
from photutils import Background2D, MedianBackground
import itertools


# function for fitting ZP equation
# this one forces slope = 1
zpfunc = lambda x, zp: x + zp

# this function allows the slope to vary
zpfuncwithslope = lambda x, m, zp: m*x + zp

pixelscale = {'HDI':0.43, 'INT':0.331, 'BOK':0.45252} 

def panstarrs_query(ra_deg, dec_deg, rad_deg, maxmag=20,
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
    field = coord.SkyCoord(ra=ra_deg, dec=dec_deg,
                           unit=(u.deg, u.deg),
                           frame='icrs')
    return vquery.query_region(field,
                               width=("%fd" % rad_deg),
                               catalog="II/349/ps1")[0]
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
    z = np.zeros_like(x)
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

def fit_circle_func():
    #fit function centered on 
    pass
class getzp():
    def __init__(self, args):
        
        self.image = args.image

        # TODO - get plot name from image header OBJECT
        header = fits.getheader(self.image)
        #self.plotprefix = self.image.split('coadd')[0].replace('.','-').replace('-noback',"")
        try:
            self.plotprefix = header['OBJECT']
        except KeyError:
            self.plotprefix = self.image.split('.fits')[0]
            
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
        self.filter = args.filter


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
        self.norder = int(args.norder)
        self.fwhm = args.fwhm


        global v1, v2
        if self.filter == 'ha':
            v1 = .9
            v2 = 1.1
        else:
            v1=.95
            v2=1.05

        
    def getzp(self):
        plt.close('all')
        if self.verbose:
            print('')
            print('STATUS: running se')        
        self.runse()
        if self.verbose:
            print('')        
            print('STATUS: getting panstarrs')
        self.get_panstarrs()
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
        
    def runse(self):
        """
        Run Source Extractor on image to measure magnitudes
        """

        os.system('ln -s ' +self.astrodir + '/default.* .')
        t = self.image.split('.fits')
        froot = t[0]
        if self.instrument == 'h':
            defaultcat = 'default.sex.HDI'
        elif self.instrument == 'i':
            defaultcat = 'default.sex.INT'
            self.keepsection=[1000,5000,0,4000]
        elif self.instrument == 'm':
            defaultcat = 'default.sex.HDI'
        elif self.instrument == 'b':
            defaultcat = 'default.sex.HDI'
        header = fits.getheader(self.image)
        try:
            expt = header['EXPTIME']
        except KeyError:
            expt = 1.
        ADUlimit = 40000.
        if self.instrument == 'i':
            if (self.filter == 'r'):
                ADUlimit = 400000./60#/float(expt)
            elif self.filter == 'ha':
                ADUlimit = 40000./180.
        #print('saturation limit in ADU/s {:.1f}'.format(ADUlimit))
        if self.fwhm is None:
            t = 'sex ' + self.image + ' -c '+defaultcat+' -CATALOG_NAME ' + froot + '.cat -MAG_ZEROPOINT 0 -SATUR_LEVEL '+str(ADUlimit)
            #t = 'sex ' + self.image + ' -c '+defaultcat+' -CATALOG_NAME ' + froot + '.cat -MAG_ZEROPOINT 0 -SATUR_LEVEL '
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
            secat_filename = froot+'.cat'
            self.secat = fits.getdata(secat_filename,2)
            self.secat0 = self.secat
            # get median fwhm of image
            # for some images, this comes back as zero, and I don't know why
            fwhm = np.median(self.secat['FWHM_IMAGE'])*self.pixelscale
            
            
            t = 'sex ' + self.image + ' -c '+defaultcat+' -CATALOG_NAME ' + froot + '.cat -MAG_ZEROPOINT 0 -SATUR_LEVEL '+str(ADUlimit)+' -SEEING_FWHM '+str(fwhm)
            if float(fwhm) == 0:
                print('WARNING: measured FWHM is zero!')
            if self.verbose:
                print('running SE again with new FWHM to get better estimate of CLASS_STAR')
        else:
            t = 'sex ' + self.image + ' -c '+defaultcat+' -CATALOG_NAME ' + froot + '.cat -MAG_ZEROPOINT 0 -SATUR_LEVEL '+str(ADUlimit)+' -SEEING_FWHM '+self.fwhm
            if self.verbose:
                print(t)
                print('running SE w/user input for FWHM to get better estimate of CLASS_STAR')            
        #############################################################
        # rerun Source Extractor catalog with updated SEEING_FWHM
        #############################################################

        #print(t)
        os.system(t)

        ###################################
        # Read in Source Extractor catalog
        ###################################

        secat_filename = froot+'.cat'
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
    def get_panstarrs(self):

        ###################################
        # get Pan-STARRS catalog over the same region
        ###################################
        ptab_name = self.image.split('.fits')[0]+'_pan_tab.csv'

        if os.path.exists(ptab_name):
            print('panstarrs table already downloaded')
            self.pan = Table.read(ptab_name)
        else:
            # check for catalog from previous run on mosaic
            print(self.image)
            header = fits.getheader(self.image)
            objname = header['OBJECT']
            filter = header['FILTER'].replace('+','').replace('nm','')
            
            glob_ptab_name = '*'+objname+'-'+filter+'_pan_tab.csv'        
            filelist = glob.glob(glob_ptab_name)
            if len(filelist) > 0:
                ptab_name = filelist[0]
                self.pan = Table.read(ptab_name)
            else:
                self.pan = panstarrs_query(self.centerRA, self.centerDEC, self.radius)
                ptab = Table(self.pan)
                ptab.write(ptab_name,format='csv',overwrite=True)
    def match_coords(self):
        ###################################
        # match Pan-STARRS1 data to Source Extractor sources
        # remove any objects that are saturated or non-linear in our r-band image
        ###################################


        secoords = SkyCoord(ra=self.secat['ALPHA_J2000'].value,dec=self.secat['DELTA_J2000'].value,unit=(u.degree,u.degree),frame='icrs')
        print(self.pan)
        pancoords = SkyCoord(ra=self.pan['RAJ2000'].value,dec=self.pan['DEJ2000'].value,unit=(u.degree,u.degree),frame='icrs')

        index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)

        # only keep matches with matched RA and Dec w/in 5 arcsec
        self.matchflag = dist2d.degree < 5./3600


        self.matchedarray1=np.zeros(len(pancoords),dtype=self.secat.dtype)
        self.matchedarray1[self.matchflag] = self.secat[index[self.matchflag]]

        ###################################################################
        # remove any objects that are saturated, have FLAGS set, galaxies,
        # must have 14 < r < 17 according to Pan-STARRS
        ###################################################################
        if self.verbose:
            print(f'\t matched {np.sum(self.matchflag)} objects')
        self.fitflag = self.matchflag  & (self.pan['rmag'] > 14.) & (self.matchedarray1['FLAGS'] <  1) & (self.pan['Qual'] < 64)  & (self.matchedarray1['CLASS_STAR'] > 0.95) & (self.pan['rmag'] < 17) #& (self.matchedarray1['MAG_AUTO'] > -11.)
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
                
        if self.filter == 'R':
            ###################################
            # Calculate Johnson R
            # from http://www.sdss3.org/dr8/algorithms/sdssUBVRITransform.php
            ###################################
            self.R = self.pan['rmag'] + (-0.153)*(self.pan['rmag']-self.pan['imag']) - 0.117

            ###################################
            # Other transformations from 
            # https://arxiv.org/pdf/1706.06147.pdf
            # R - r = C0 + C1 x (r-i)  (-0.166, -0.275)
            # R - r = C0 + C1 x (g-r)  (-0.142, -0.166)
            ###################################
            #
            if self.useri:
                self.R = self.pan['rmag'] + (-0.166)*(self.pan['rmag']-self.pan['imag']) - 0.275
            else:
                self.R = self.pan['rmag'] + (-0.142)*(self.pan['gmag']-self.pan['rmag']) - 0.142

        else:
            self.R = self.pan['rmag']
            
    def plot_fitresults(self, x, y, yerr=None, polyfit_results = [0,0]):
        # plot best-fit results
        yfit = np.polyval(polyfit_results,x)
        residual = (yfit - y)
        plt.figure(figsize=(8,8))
        s = ' (MAD = %.2f)'%(MAD2(residual))
        plt.title(self.plotprefix+s)
        
        plt.subplot(2,1,1)
        if len(yerr) < len(y): 
            plt.plot(x,y,'bo',label='MAG_AUTO')
            
        else:
            plt.errorbar(x,y,yerr=yerr,fmt='bo',ecolor='b',label='SE MAG')
        plt.xlabel('Pan-STARRS r',fontsize=16)
        plt.ylabel('SE R-band MAG',fontsize=16)
        xl = np.linspace(14,17,10)
        yl = np.polyval(polyfit_results,xl)
        s = 'fit: y = %.2f PAN + %.2f'%(polyfit_results[0],polyfit_results[1])
        plt.plot(xl,yl,'k--',label=s)
        plt.legend()
        
        plt.subplot(2,1,2)
        s = 'MAD = %.4f'%(MAD2(residual))
        if len(yerr) < len(y):
            plt.plot(x,residual, 'ko',label=s)
            
        else:
            plt.errorbar(x,residual,yerr=yerr,fmt='None',ecolor='b',label='SE MAG '+s)
        plt.xlabel('Pan-STARRS r',fontsize=16)
        plt.ylabel('YFIT - SE R-band MAG',fontsize=16)
        plt.legend()
        plt.axhline(y=0,color='r')
        plt.savefig('plots/'+self.plotprefix.replace('.fits','')+'se-pan-flux.png')
        plt.close()
    def fitzp(self,plotall=False):
        ###################################
        # Solve for the zeropoint
        ###################################

        # plot Pan-STARRS r mag on x axis, observed R-mag on y axis
        flag = self.fitflag
        c = np.polyfit(self.pan['rmag'][flag],self.matchedarray1['MAG_AUTO'][flag],1)

        if plotall:
            plt.figure(figsize=(6,4))
            plt.title(self.plotprefix)
            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_AUTO'][flag],'bo')
            plt.errorbar(self.pan['rmag'][flag],self.matchedarray1['MAG_AUTO'][flag],xerr= self.pan['e_rmag'][flag],yerr=self.matchedarray1['MAGERR_AUTO'][flag],fmt='none')
            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_BEST'][flag],'ro',label='MAG_BEST')
            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_PETRO'][flag],'go',label='MAG_PETRO')
            plt.plot(self.pan['rmag'][flag],self.matchedarray1['MAG_APER'][:,0][flag],'ko',label='MAG_APER')
            plt.xlabel('Pan-STARRS r',fontsize=16)
            plt.ylabel('SE R-band MAG_AUTO',fontsize=16)

            xl = np.linspace(14,17,10)
            yl = np.polyval(c,xl)
            plt.plot(xl,yl,'k--')
            #plt.savefig('getzp-fig2.png')
            #plt.plot(xl,1.2*yl,'k:')
            #print(c)
    
        yfit = np.polyval(c,self.pan['rmag'])
        residual = np.zeros(len(flag))
        ####################################
        ## had been dividing by yfit, but that doesn't make sense
        ## want residual to be in magnitudes
        ## removing yfit normalization
        ####################################
        residual[flag] = (yfit[flag] - self.matchedarray1['MAG_AUTO'][flag])#/yfit[flag]

        self.bestc = np.array([0,0],'f')
        delta = 100.     
        x = self.R[flag] # expected mag from panstarrs
        # fixed radii apertures: [:,0] = 3 pix, [:,1] = 5 pix, [:,2] = 7 pixels

        if self.mag == 0: # this is the default magnitude
            if self.verbose:
                print('Using Aperture Magnitudes')
            y = self.matchedarray1['MAG_APER'][:,self.naper][flag]
            yerr = self.matchedarray1['MAGERR_APER'][:,self.naper][flag]
        elif self.mag == 1:
            if self.verbose:
                print('Using MAG_BEST')
            y = self.matchedarray1['MAG_BEST'][flag]
            yerr = self.matchedarray1['MAGERR_BEST'][flag]
        elif self.mag == 2:
            if self.verbose:
                print('Using MAG_PETRO')
            y = self.matchedarray1['MAG_PETRO'][flag]
            yerr = self.matchedarray1['MAGERR_PETRO'][flag]
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
                self.plot_fitresults(x,y,yerr=yerr,polyfit_results = c)

    
            # check for convergence
            if self.verbose:
                print('new ZP = {:.3f}, previous ZP = {:.3f}'.format(self.bestc[1],c[1]))
            delta = abs(self.bestc[1] - c[1])
            self.bestc = c
            MAD = 1.48*np.median(abs(residual - np.median(residual)))
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
        ###################################
        ##  show histogram of residuals
        ###################################

        yplot = self.matchedarray1['MAG_APER'][:,self.naper][self.fitflag]
        magfit = np.polyval(self.bestc,self.R[self.fitflag])
        residual_all = 10.**((magfit - yplot)/2.5)        
        self.residual_all = residual_all
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
        
        s = ' (mean,std,MAD = {:.2f},{:.2f},{:.2f})'.format(np.mean(residual_all),np.std(residual_all),MAD2(residual_all))
        #s = str(MAD(residual_all))
        plt.title(self.plotprefix+s)
        self.residual_allx = self.matchedarray1['X_IMAGE'][self.fitflag]
        self.residual_ally = self.matchedarray1['Y_IMAGE'][self.fitflag]
        plt.scatter(self.matchedarray1['X_IMAGE'][self.fitflag],self.matchedarray1['Y_IMAGE'][self.fitflag],c = (residual_all),vmin=v1,vmax=v2,s=15)
        cb=plt.colorbar()
        cb.set_label('f-meas/f-pan')
        plt.savefig('plots/'+self.plotprefix.replace(".fits","")+'getzp-xyresidual-fitted.png')

        self.x = x
        self.y = y
        self.yerr = yerr
        self.zpcovar = t[1]
        self.zperr = np.sqrt(self.zpcovar[0][0])
        self.zp = self.bestc[1]
        self.plot_fitresults(x,y,yerr=yerr,polyfit_results = self.bestc)
                
    def update_header(self):
        #print('working on this')
        # add best-fit ZP to image header
        im, header = fits.getdata(self.image,header=True)
        zperr = np.sqrt(self.zpcovar[0][0])            
        # or convert vega zp to AB
        if self.filter == 'R':
            # conversion from Blanton+2007
            # http://www.astronomy.ohio-state.edu/~martini/usefuldata.html
            header.set('PHOTZP',float('{:.3f}'.format(-1.*self.bestc[1]+.21)))

            header.set('LAMB(um)',float(.6442))

        else:
            header.set('PHOTZP',float('{:.3f}'.format(-1.*self.bestc[1])))
        header.set('PHOTZPER',float('{:.3f}'.format(zperr)))            
        header.set('PHOTSYS','AB')
        header.set('FLUXZPJY',float(3631))
        fits.writeto(self.image, im, header, overwrite=True)
        
    def fit_residual_surface(self,norder=2,suffix=None):
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
        self.weightdata = fits.getdata(self.image)
        # this is not used
        self.nodata =  self.weightdata == 0

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
        
        # Fit a 2nd order, 2d polynomial
        m = polyfit2d(self.xim[~clip_flag.mask],self.yim[~clip_flag.mask],self.zim[~clip_flag.mask],order=norder)

        # Evaluate it on a grid...        
        ny, nx = self.imagedata.shape
        xx, yy = np.meshgrid(np.linspace(self.xim.min(), self.xim.max(), nx), 
                         np.linspace(self.yim.min(), self.yim.max(), ny))
        self.zz = polyval2d(xx, yy, m)

        # Plot
        plt.figure()
        plt.imshow(self.zz,extent=(self.xim.min(), self.yim.max(), self.xim.max(), self.yim.min()),vmin=v1,vmax=v2)
        plt.scatter(self.xim[~clip_flag.mask], self.yim[~clip_flag.mask], c=self.zim[~clip_flag.mask],vmin=v1,vmax=v2,s=15)
        cb=plt.colorbar()
        cb.set_label('f-meas/f-pan')
        s = ' std (MAD) = %.4f (%.4f)'%(np.std(self.zim[~clip_flag.mask]),MAD2(self.zim[~clip_flag.mask]))
        plt.title(self.plotprefix+': n poly = '+str(norder)+s)
        #plt.show()
        
        if suffix is None:
            plotname='-imsurfit-'+str(norder)
        else:
            plotname='-imsurfit-'+str(norder)+'-'+suffix
        plt.savefig('plots/'+self.plotprefix.replace(".fits","")+plotname+'.png')
        plt.savefig('plots/'+self.plotprefix.replace(".fits","")+plotname+'.pdf')
        
    def renorm_wfc(self):
        # normalize surface fit
        self.zz_norm = self.zz/np.median(self.zz)
        # not sure we want to normalize this, actually
        self.zz_norm = self.zz
        
        # divide image by surface fit
        self.imagedata,header = fits.getdata(self.image,header=True)
        self.imagedata_norm = self.imagedata/self.zz_norm
        # save flattened image
        self.renorm_image = 'f'+self.image
        fits.writeto(self.renorm_image,self.imagedata_norm,header=header,overwrite=True)
    def rerun_zp_fit(self):
        # change image name to flattened image
        self.image = self.renorm_image
        self.plotprefix = 'f'+self.plotprefix
        # rerun getzp, but don't download panstarrs again
        self.runse()
        if self.verbose:
            print('STATUS: matching se cat to panstarrs')       
        self.match_coords()
        if self.verbose:
            print('STATUS: fitting zeropoint')        
        self.fitzp()
        if self.verbose:
            print('STATUS: udating header')        
        self.update_header()
        # check to make sure the systematic residuals have been removed
        self.fit_residual_surface(suffix='round2',norder=self.norder)


def main(raw_args=None):
    parser = argparse.ArgumentParser(description ='Run sextractor, get Pan-STARRS catalog, and then computer photometric ZP\n \n from within ipython: \n %run ~/github/Virgo/programs/getzp.py --image pointing031-r.coadd.fits --instrument i \n \n The y intercept is -1*ZP. \n \n x and y data can be accessed at zp.x and zp.y in case you want to make additional plots.', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('--image', dest = 'image', default = 'test.coadd.fits', help = 'Image for ZP calibration')
    parser.add_argument('--instrument', dest = 'instrument', default = None, help = 'HDI = h, KPNO mosaic = m, INT = i, BOK 90Prime = b')
    parser.add_argument('--catalog', dest = 'catalog', default = None, help = 'photometric catalog to use for bootrapping photometry')    
    parser.add_argument('--fwhm', dest = 'fwhm', default = None, help = 'image FWHM in arcseconds.  Default is none, then SE assumes 1.5 arcsec')    
    parser.add_argument('--filter', dest = 'filter', default = 'R', help = 'filter (R or r; use ha for Halpha)')
    parser.add_argument('--useri',dest = 'useri', default = False, action = 'store_true', help = 'Use r->R transformation as a function of r-i rather than the g-r relation.  g-r is the default.')
    parser.add_argument('--normbyexptime', dest = 'normbyexptime', default = False, action = 'store_true', help = "set this flag if the image is in ADU rather than ADU/s, and the program will then normalize by the exposure time.  Note: swarp produces images in ADU/s, so this is usually not necessary if using coadds from swarp.")
    parser.add_argument('--mag', dest = 'mag', default = 0,help = "select SE magnitude to use when solving for ZP.  0=MAG_APER,1=MAG_BEST,2=MAG_PETRO.  Default is MAG_APER ",choices=['0','1','2'])
    parser.add_argument('--naper', dest = 'naper', default = 5,help = "select fixed aperture magnitude.  0=10pix,1=12pix,2=15pix,3=20pix,4=25pix,5=30pix.  Default is 5 (30 pixel diameter)")
    parser.add_argument('--nsigma', dest = 'nsigma', default = 3.5, help = 'number of std to use in iterative rejection of ZP fitting.  default is 3.5')
    parser.add_argument('--d',dest = 'd', default ='~/github/HalphaImaging/astromatic/', help = 'Locates path of default config files.  Default is ~/github/HalphaImaging/astromatic')
    parser.add_argument('--fit',dest = 'fitonly', default = False, action = 'store_true',help = 'Do not run SE or download catalog.  just redo fitting.')
    parser.add_argument('--flatten',dest = 'flatten', default = 0, help = 'Number of time to run flattening process to try to remove vignetting/illumination patterns.  The default is zero.  Options are [0,1,2].  This is needed for INT data from 2019.  HDI does not show this effect, and INT data from 2022 does not seem to show it either.',choices=['0','1','2'])    
    parser.add_argument('--norder',dest = 'norder', default = 2, help = 'degree of polynomial to fit to overall background.  default is 2.',choices=['0','1','2'])    
    parser.add_argument('--verbose',dest = 'verbose', default = False, action = 'store_true',help = 'print extra debug/status statements')
    parser.add_argument('--getrefcatonly',dest = 'getrefcatonly',default=False,action='store_true',help='download reference PANSTARRS catalog only.  use this before running with slurm')
    args = parser.parse_args(raw_args)
    #zp = getzp(args.image, instrument=args.instrument, filter=args.filter, astromatic_dir = args.d,norm_exptime = args.nexptime, nsigma = float(args.nsigma), useri = args.useri,naper = args.naper, mag = int(args.mag))
    zp = getzp(args)
    zp.getzp()
    print('ZP = {:.3f} +/- {:.3f}, {}'.format(-1*zp.zp,zp.zperr,zp.image))
    return zp,-1*zp.zp,zp.zperr


if __name__ == '__main__':
    main()
