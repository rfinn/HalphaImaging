#!/usr/bin/env python

"""
GOAL:
* will download images from legacy server, unwise, and galex
* plot cutouts together

INPUT:
* image name of r-band cutout
-OR-
* RA, DEC, size and ID

OUTPUT:
can create different version of cutouts - halpha, sfr indicators, and all

USAGE:
* Images will get saved to the current directory, so run this from individual cutout directories.
* run this after the halpha cutouts have been made
* see main below

UPDATES ON 7/7/2020:
* trying to adapt this to take an RA, DEC and size instead of the r-band cutout
  - goal is to get cutouts for ALL VF galaxies

STILL TO DO:
* need to stack unwise images when multiple images are returned

REFERENCES
* For stacking unWISE images:
https://reproject.readthedocs.io/en/stable/celestial.html

"""


import sys
import wget
import tarfile
import os
import numpy as np
import glob
import argparse
from argparse import RawDescriptionHelpFormatter

from matplotlib import pyplot as plt
from PIL import Image
from scipy.stats import scoreatpercentile

from urllib.parse import urlencode
from urllib.request import urlretrieve

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord

from astropy import units as u
from astropy.nddata import Cutout2D
from astropy.stats import sigma_clip
from astropy.visualization import simple_norm
from astroquery.mast import Observations

vmin = .5
vmax = 50.

UNWISE_PIXSCALE = 2.75
LEGACY_PIXSCALE = 1


def get_legacy_images(ra,dec,galid='VFID0',pixscale=1,imsize='60',band='g',makeplots=False,subfolder=None):
    """
    Download legacy image for a particular ra, dec
    
    Inputs:
    * ra
    * dec
    * galid = galaxy id (e.g. VFID0001); used for naming the image files
    * imsize = size of cutout in pixels
    * band = filter for the fits images that will be returned.  e.g. 'g' or 'r' or 'z'. 
    * pixscale = pixel scale of cutout in arcsec; native is 0.262 for legacy
    * makeplots = boolean, generate plot of image
    * subfolder = default is None; you can specify a name of a subfolder to 
                  save the data in, e.g., subfolder='legacy-images'
    Returns:
    * fits_name = fits image name
    * jpeg_name = jpeg image name
    """
    imsize = int(imsize)

    # make output image names
    if subfolder is not None:
        # check if subfolder exists. if not, make it.
        if not os.path.exists(subfolder):
            os.mkdir(subfolder)
        rootname = subfolder+'/'+str(galid)+'-legacy-'+str(imsize)
    else:
        rootname = str(galid)+'-legacy-'+str(imsize)        
    jpeg_name = rootname+'.jpg'
    fits_name = rootname+'-'+band+'.fits'


    print('legacy imsize = ',imsize)
    
    # check if images already exist
    # if not download images
    if not(os.path.exists(jpeg_name)):
        print('retrieving ',jpeg_name)
        url='http://legacysurvey.org/viewer/jpeg-cutout?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(pixscale)
        urlretrieve(url, jpeg_name)
    else:
        print('previously downloaded ',jpeg_name)
    if not(os.path.exists(fits_name)):
        print('retrieving ',fits_name)
        url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra)+'&dec='+str(dec)+'&layer=dr8&size='+str(imsize)+'&pixscale='+str(pixscale)+'&bands='+band
        urlretrieve(url, fits_name)
    else:
        print('previously downloaded ',fits_name)

    # try to read the data in
    try:
        t,h = fits.getdata(fits_name,header=True)
        
    except IndexError:
        print('problem accessing image')
        print(fits_name)
        url='http://legacysurvey.org/viewer/cutout.fits?ra='+str(ra1)+'&dec='+str(dec1)+'&layer=dr8&size='+str(image_size)+'&pixscale=1.00'
        print(url)
        return None

    # this will trigger if the image is outside the legacy footprint
    if np.mean(t[1]) == 0:
        return None

    # plot the images
    if makeplots:
        if jpeg:
            t = Image.open(jpeg_name)
            plt.imshow(t,origin='upper')
        else:
            norm = simple_norm(t[1],stretch='asinh',percent=99.5)            
            plt.imshow(t[1],origin='upper',cmap='gray_r', norm=norm)

    # return the name of the fits images and jpeg image
    return fits_name, jpeg_name


def get_unwise_image(ra,dec,galid='VFID0',pixscale=2.75,imsize='60',bands='1234',makeplots=False,subfolder=None):
    """
    Download unwise image for a particular ra, dec
    
    Inputs:
    * ra
    * dec
    * galid = galaxy id (e.g. VFID0001)
    * imsize = size of cutout in pixels
    * pixscale = pixel scale of cutout in arcsec
      - native is 0.262 for legacy; 
      - 2.75 for wise
    """
    downloadwise = True
    # check if images already exist
    if subfolder is not None:
        image_names = glob.glob(subfolder+'/'+galid+'-unwise*img-m.fits')
    else:
        image_names = glob.glob(galid+'-unwise*img-m.fits')
    if len(image_names) > 3:
        print('unwise images already downloaded')
        print(image_names)
        # should be only one *-img-m.fits image per band
        if len(image_names) > len(bands):
            multiframe=True
        else:
            multiframe = False
        weight_names = glob.glob(galid+'-unwise*std*.fits')
        if not multiframe:
            return image_names,weight_names,multiframe
        else:
            print('going to try new stacking for wise multiframe')
            downloadwise = False
    if downloadwise:
        imsize = int(imsize)
        print('wise image size = ',imsize)
        baseurl = 'http://unwise.me/cutout_fits?version=allwise'
        imurl = baseurl +'&ra=%.5f&dec=%.5f&size=%s&bands=%s'%(ra,dec,imsize,bands)
        print('downloading unwise images')
        print(imurl)
        wisetar = wget.download(imurl)
        tartemp = tarfile.open(wisetar,mode='r:gz') #mode='r:gz'
        wnames = tartemp.getnames()

        print(wnames)
        # check for multiple pointings - means galaxy is split between images
        multiframe = False
        if len(wnames) > 4*len(bands):
            multiframe = True
    
        wmembers = tartemp.getmembers()
        image_names = []
        weight_names = []
        tartemp.extractall()
        for fname in wnames:
            #print(fname)
            t = fname.split('-')
            if subfolder is not None:
                rename = subfolder+'/'+str(galid)+'-'+t[0]+'-'+t[1]+'-'+t[2]+'-'+t[3]+'-'+t[4]
            else:
                rename = str(galid)+'-'+t[0]+'-'+t[1]+'-'+t[2]+'-'+t[3]+'-'+t[4]
            #print('rename = ',rename)
            #print(rename.find('gz'))
            #if os.path.exists(rename): # this should only occur if multiple images are returned from wise
            #    os.remove(rename)
            os.rename(fname, rename)
            if rename.find('.gz') > -1:
                #print('hello????')
                os.system('gunzip '+rename)
                rename = rename.split('.gz')[0]
                #print('after gunzip, rename = ',rename)
            if rename.find('img') > -1:
                image_names.append(rename)
            if rename.find('std') > -1:
                # move ivar images to imagename.weight.fits
                outim = rename.replace('std-m','img-m.std')
                os.rename(rename,outim)
                weight_names.append(outim)
        os.remove(wisetar)
    
    # if multiframe
    # run swarp to create coadded image
    if multiframe:
        image_names=[]
        weight_names=[]
        for b in bands:
            print('running swarp to combine multiple unwise images in band ',b)
            #########################################
            ## COMBINE THE IMAGE FRAMES USING AVERAGE
            #########################################        
            # create default.swarp
            os.system('swarp -d > default.swarp')
            # run swarp
            matchstring = "*w{}-img-m.fits".format(b)            
            if subfolder is not None:
                allfiles = glob.glob(subfolder+'/'+galid+matchstring)


                
            else:
                allfiles = glob.glob(galid+matchstring)

            all_images = " ".join(allfiles)
            output_image = str(galid)+'-'
            s = 'swarp '+all_images+' -COMBINE_TYPE AVERAGE -WEIGHT_SUFFIX .std.fits -SUBTRACT_BACK N'
            print(s)
            os.system(s)

            # rename coadd.fits to the output image name
            outimage = str(galid)+'-unwise-w'+str(b)+'-coadd.fits'
            if subfolder is not None:
                os.rename('coadd.fits',os.path.join(subfolder,outimage))
                image_names.append(os.path.join(subfolder,outimage))                
            else:
                os.rename('coadd.fits',outimage)
                image_names.append(outimage)

            
            #os.rename('coadd.fits','unwise/'+outimage)

            #########################################
            ## COMBINE THE STD FRAMES USING SUM
            #########################################
            # just combine using average and then multiply by sqrt 2
            matchstring = "*w{}-img-m.std.fits".format(b)
            
            if subfolder is not None:
                allfiles = glob.glob(subfolder+'/'+galid+matchstring)
            else:
                allfiles = glob.glob(galid+matchstring)
            #print('allfiles with std images = ',allfiles)
            all_images = " ".join(allfiles)
            
            output_image = str(galid)+'-'
            s = 'swarp '+all_images+' -COMBINE_TYPE MEDIAN -WEIGHT_TYPE NONE -SUBTRACT_BACK N'
            #print(s)
            os.system(s)

            # now take sqrt of image values
            #im,h = fits.getdata('coadd.fits',header=True)
            #im = np.sqrt(2)*im
            # divide by number of images because we took average of data
            #im = im/len(all_images)
            
            weightname = str(galid)+'-unwise-w'+str(b)+'-coadd.std.fits'  
            #fits.writeto(weightname,im,header=h,overwrite=True)


            # just trying average!
            if subfolder is not None:
                os.rename('coadd.fits',os.path.join(subfolder,weightname))
                weight_names.append(os.path.join(subfolder,weightname))
            else:
                os.rename('coadd.fits',weightname)
                weight_names.append(weightname)
            
        
    if makeplots:
        ##### DISPLAY IMAGE
        im = fits.getdata(rename)
        norm = simple_norm(im, stretch='asinh',percent=99)
        plt.imshow(im, norm=norm,origin='upper')
        #plt.show()
    #print(image_names)
    #print(multiframe)


    return image_names,weight_names,multiframe

def get_galex_image(ra,dec,imsize):
    """

    get galex image of a galaxy
    
    Input:
    * ra in deg
    * dec in deg
    * imsize in arcsec
    
    Returns:
    * image
    """
    

    # following procedure outlined here:
    # https://astroquery.readthedocs.io/en/latest/mast/mast.html

    # get data products in region near ra,dec
    obs_table = Observations.query_region("%12.8f %12.8f"%(ra,dec),radius=.1*u.arcmin)
    # create a flag to select galex data
    galexFlag = obs_table['obs_collection'] == 'GALEX'

    # separate out galex data
    data_products = Observations.get_product_list(obs_table[galexFlag])

    # download the observations
    manifest = Observations.download_products(data_products,productType="SCIENCE")

    for m in manifest:
        # choose the first NUV image
        if m['Local Path'].find('nd-int') > -1:
            nuv_path = m['Local Path']
            break
        
    # should be able to construct path from the obs_id in data_products
    # this will let us check if the image is already downloaded
    #
    # DONE: I can also save the cutout in a GALEX folder, and
    # look for the image before calling this function
    #
    # but I should also look for the image, because if I change the image
    # size of the cutout, I don't need to download the big FOV
    
    nuv,nuv_header = fits.getdata(nuv_path,header=True)

    # this is a big image, so we need to get a cutout

    nuv_wcs = WCS(nuv_header)
    position = SkyCoord(ra,dec,unit="deg",frame='icrs')
    print("\nYoohoo!\n")
    cutout = Cutout2D(nuv,position,(imsize*u.arcsec,imsize*u.arcsec),wcs=nuv_wcs)
    #try:
    #    cutout = Cutout2D(nuv,position,(imsize*u.arcsec,imsize*u.arcsec),wcs=nuv_wcs)
    #except:
    #    print('WARNING: problem getting galex cutout')
    #    cutout = None

    return cutout
    
def display_image(image,percent=99.9,lowrange=False,mask=None,sigclip=True):
    lowrange=False
    if sigclip:
        clipped_data = sigma_clip(image,sigma_lower=5,sigma_upper=5)#,grow=10)
    else:
        clipped_data = image
    if lowrange:
        norm = simple_norm(clipped_data, stretch='linear',percent=percent)
    else:
        norm = simple_norm(clipped_data, stretch='asinh',percent=percent)

    plt.imshow(image, norm=norm,cmap='gray_r',origin='lower')
    #v1,v2=scoreatpercentile(image,[.5,99.5])            
    #plt.imshow(image, cmap='gray_r',vmin=v1,vmax=v2,origin='lower')    

class cutouts():
    def __init__(self, rimage,ra=None,dec=None,size=None,galid=None):
        '''
        arguments:
        - rimage - rband cutout image; position and image size will be taken from this
        - ra - alternatively, user can provide ra, dec, size, and rootname
        - dec - dec of center of galaxy or center of cutout
        - size - size of cutout in arcsec
        - galid - prefix to use with all the cutouts; should be VFIDXXXX-nedname-
        '''
        self.galid = galid
        if ra is not None:
            self.ra = ra
            self.dec = dec
            self.xsize_arcsec = size
            self.ysize_arcsec = size            
            self.galid = galid
            self.rootname = galid
            self.use_position_flag = True
            # look for r-band cutout
            self.find_rcutout()
        else:
            self.use_position_flag = False
            self.r_name = rimage
            self.wcs = WCS(self.r_name)

            if self.r_name.find('_R') > -1:
                split_string = '_R.fits'

            if self.r_name.find('-R') > -1:
                split_string = '-R.fits'
            elif self.r_name.find('-r') > -1:
                split_string = '-r.fits'
            elif self.r_name.find('_r') > -1:
                split_string = '_r.fits'
            self.rootname = self.r_name.split(split_string)[0]
        # check to make sure cutouts directory exists
        if not os.path.exists('galex'):
            os.mkdir('galex')
        if not os.path.exists('legacy'):
            os.mkdir('legacy')
        if not os.path.exists('unwise'):
            os.mkdir('unwise')
                        
    def runall(self):
        '''
        get halpha cutout if it exists

        if using r-band image to as reference, get image size, ra, dec, id
        '''
        try:
            self.get_halpha_cutouts()
            self.halpha_flag = True
        except:
            print(self.galid,': did not find halpha cutout')
            self.halpha_flag = False            
        if self.use_position_flag == False:
            self.get_image_size()
            self.get_RADEC()
            self.get_galid()


        self.get_galex_image()
        try:
            self.get_galex_image()
        except:
            print("\nWARNING: Problem with GALEX query - maybe outside footprint?\n")
            self.nuv_flag = False
            self.galex_image = None
        self.download_unwise_images()
        self.load_unwise_images()

        try:
            print('trying to download legacy image')
            self.download_legacy()
            self.load_legacy_images()
            self.legacy_flag = True
        except: # urllib.error.HTTPError:
            print('\nWARNING: could not get legacy image\n')
            self.legacy_flag = False
        

        #self.plotallcutouts()
    def runha(self):
        self.get_halpha_cutouts()
        self.get_image_size()
        self.get_RADEC()
        self.get_galid()
        self.plotcutouts()
    def find_rcutout(self):
        imlist = glob.glob(self.galid+'*-R.fits')
        if len(imlist) == 0:
            imlist = glob.glob(self.galid+'*-r.fits')
        if len(imlist) == 0:
            self.ha_flag = False
            self.r_name = None
        else:
            self.ha_flag = True
            self.r_name = imlist[0]
        
    def get_halpha_cutouts(self):
        if self.r_name is not None:
            self.r,self.header = fits.getdata(self.r_name,header=True)
            self.ha = fits.getdata(self.rootname+'-Ha.fits')
            self.cs = fits.getdata(self.rootname+'-CS.fits')
            try:
                self.mask = fits.getdata(self.rootname+'-R-mask.fits')
            except:
                print('WARNING: could not open mask ',self.rootname+'-R-mask.fits')
    def get_image_size(self):
        # get image size in pixels and arcsec
        self.xsize_pix,self.ysize_pix = self.r.shape
        try:
            self.pscale = np.abs(float(self.header['PIXSCAL1']))/3600
        except KeyError:
            self.pscale = np.abs(float(self.header['CD1_1']))
        self.xsize_arcsec = self.xsize_pix*self.pscale*3600
        self.ysize_arcsec = self.ysize_pix*self.pscale*3600    
    def get_RADEC(self):
        xcenter = self.xsize_pix/2
        ycenter = self.ysize_pix/2
        self.ra,self.dec = self.wcs.wcs_pix2world(xcenter,ycenter,1)
    def get_galid(self):
        try:
            self.galid = self.header['ID']
        except:
            self.galid = self.galid
        
    def download_legacy(self):

        # get legacy grz color and fits
        self.legacy_imsize = self.xsize_arcsec/LEGACY_PIXSCALE
        print('requested legacy imsize = ',self.legacy_imsize)
        self.legacy_filename_g,self.legacy_jpegname = get_legacy_images(self.ra,self.dec,galid=self.galid,imsize=self.legacy_imsize,band='g',subfolder='legacy')
        h = fits.getheader(self.legacy_filename_g)
        self.legacy_pscale = np.abs(h['CD1_1'])*3600
        self.legacy_filename_r,self.legacy_jpegname = get_legacy_images(self.ra,self.dec,galid=self.galid,imsize=self.legacy_imsize,band='r',subfolder='legacy')
        self.legacy_filename_z,self.legacy_jpegname = get_legacy_images(self.ra,self.dec,galid=self.galid,imsize=self.legacy_imsize,band='z',subfolder='legacy')
        
    def load_legacy_images(self):
        self.legacy_g = fits.getdata(self.legacy_filename_g)
        self.legacy_r = fits.getdata(self.legacy_filename_r)
        self.legacy_z = fits.getdata(self.legacy_filename_z)        
    def download_unwise_images(self,band='1234'):
        '''
        GOAL: Get the unWISE image from the unWISE catalog

        INPUT: ra,dec used to grab unwise image information

        OUTPUT: Name of file to retrieve from

        '''

        self.wise_band = band
        self.wise_imsize = self.xsize_arcsec/UNWISE_PIXSCALE
        print('wise image size = ',self.wise_imsize)
        # check for coadds in unwise folder.  if they are there, then use these
        outimage = str(self.galid)+'-unwise-w*-coadd.fits'
        coadds = glob.glob('unwise/'+outimage)
        if len(coadds) > 0:
            print('found unwise coadds')
            self.wise_filenames = [os.path.basename(f) for f in coadds]
            stdimage = str(self.galid)+'-unwise-w*-coadd.std.fits'
            t = glob.glob('unwise/'+stdimage)
            self.wise_weightimages = [os.path.basename(f) for f in t]
            self.wise_multiframe_flag = True
        else:
            self.wise_filenames,self.wise_weightimages,self.wise_multiframe_flag = \
                get_unwise_image(self.ra,self.dec,galid=self.galid,pixscale='2.75',imsize=self.wise_imsize,bands=self.wise_band,subfolder='unwise')
            
    def load_unwise_images(self):
        print('current directory = ',os.getcwd())                
        if self.wise_multiframe_flag:
            print('WARNING: galaxy falls on multiple unwise images')
        print(self.wise_filenames)
        for f in self.wise_filenames:
            #f = os.path.join('unwise',f)

            # this is a kludge to get around issue that without multiframe
            # the file has the path prepended
            # but with multiframe, it doesn't
            #
            # but then when don't I just update filename in get_unwise_image???
            if not os.path.exists(f):
                f = 'unwise/'+f
            print('wise filename :',f)            
            if f.find('w1') > -1:
                self.w1,self.w1_header = fits.getdata(f,header=True)
                self.w1_pscale = np.abs(self.w1_header['CD1_1'])*3600
            elif f.find('w2') > -1:
                self.w2,self.w2_header = fits.getdata(f,header=True)
                self.w2_pscale = np.abs(self.w2_header['CD1_1'])*3600
            elif f.find('w3') > -1:
                self.w3,self.w3_header = fits.getdata(f,header=True)
                self.w3_pscale = np.abs(self.w3_header['CD1_1'])*3600
            elif f.find('w4') > -1:
                self.w4,self.w4_header = fits.getdata(f,header=True)
                self.w4_pscale = np.abs(self.w4_header['CD1_1'])*3600
    def get_galex_image(self):
        # download galex image if necessary (large FOV)
        # then get cutout

        imsize_arcsec = "%i"%(self.xsize_arcsec)        
        try:
            t = self.galid.split('-')
            self.nuv_image_name = 'galex/'+self.galid+'-'+t[1]+'-nuv-'+imsize_arcsec+'.fits'            
        except IndexError:
            self.nuv_image_name = 'galex/'+self.galid+'-nuv-'+imsize_arcsec+'.fits'            


        
        if os.path.exists(self.nuv_image_name):
            self.nuv_image,h = fits.getdata(self.nuv_image_name,header=True)
            self.nuv_flag = True
            try:
                self.nuv_pscale = np.abs(h['CD1_1'])*3600
            except:
                print('WARNING: could not get galex pixelscale')
        else:
            cutout = get_galex_image(self.ra,self.dec,self.xsize_arcsec)
            if cutout is not None:
                
                fits.writeto(self.nuv_image_name, cutout.data, header=cutout.header, overwrite=True)
                self.nuv_image = cutout.data
                self.nuv_flag = True
            else:
                self.nuv_flag = False
    def plotcutouts(self,plotsingle=True):
        if plotsingle:
            figure_size=(10,4)
            plt.figure(figsize=figure_size)
            plt.clf()
            plt.subplots_adjust(hspace=0,wspace=0)
        plt.subplot(1,3,1)
        self.plot_ha()
        plt.subplot(1,3,2)
        self.plot_r()
        plt.subplot(1,3,3)
        self.plot_cs()
        #plt.show(block=False)
        if plotsingle:
            plt.savefig(self.rootname+'-cutouts.png')
            plt.savefig(self.rootname+'-cutouts.pdf')
    def plotmask(self,plotsingle=True):
        if plotsingle:
            figure_size=(10,4)
            plt.figure(figsize=figure_size)
            plt.clf()
            plt.subplots_adjust(hspace=0,wspace=0)
        plt.subplot(1,3,1)
        self.plot_r()
        plt.subplot(1,3,2)
        self.plot_cs()
        plt.subplot(1,3,3)
        self.plot_mask()
        #plt.show(block=False)        
        plt.savefig(self.rootname+'-mask.png')
        plt.savefig(self.rootname+'-mask.pdf')
        
    def plotallcutouts(self,plotsingle=True):
        '''plot all wavelengths'''
        nrow = 3
        ncol = 4

        #Continuum subtracted image
        

        if plotsingle:
            figure_size=(9,7)
            plt.figure(figsize=figure_size)
            plt.clf()
            plt.subplots_adjust(left=.05,right=.95,bottom=.05,top=.9,hspace=.275,wspace=0)
        for i in range(nrow*ncol):
            
            plt.subplot(nrow,ncol,i+1)
            if i == 0:
                if self.legacy_flag:
                    self.plot_legacy_jpg()
            elif i < 4:
                # plot each band of legacy
                if self.legacy_flag:
                    self.plot_legacy(band=i)
            elif (i > 3) & (i < 8):
                wband = i-3
                self.plot_unwise(band=wband)
            elif i == 8:
                try:
                    self.plot_ha()
                except:
                    print('no ha')
            elif i == 9:
                try:
                    self.plot_r()
                except:
                    print('no r')
            elif i == 10:
                try:
                    self.plot_cs()
                except:
                    print('no cs halpha')
            elif i == 11:
                if self.nuv_flag:
                    self.plot_galex_nuv()
            #if (i%4 != 0):
            if (i > 7) & (i < 11):
                ax = plt.gca()
                mycolor = 'r'
                ax.spines['bottom'].set_color(mycolor)
                ax.spines['top'].set_color(mycolor) 
                ax.spines['right'].set_color(mycolor)
                ax.spines['left'].set_color(mycolor)
            #plt.gca().set_yticks(())
            #plt.gca().set_xticks(())            
        plt.text(-1.3,3.8,self.galid,transform=plt.gca().transAxes,fontsize=14,horizontalalignment='center')    
        plt.savefig(self.galid+'-all-cutouts.png')
        plt.savefig(self.galid+'-all-cutouts.pdf')        

    def plotsfrcutouts(self,plotsingle=True):
        nrow = 1
        ncol = 4

        #Continuum subtracted image
        

        if plotsingle:
            figure_size=(12,4)
            plt.figure(figsize=figure_size)
            plt.clf()
            plt.subplots_adjust(left=.05,right=.95,bottom=.05,top=.9,hspace=.2,wspace=0)
        for i in range(nrow*ncol):
            
            plt.subplot(nrow,ncol,i+1)
            if i == 0:
                if self.legacy_flag:
                    self.plot_legacy_jpg()
            elif i == 1:
                self.plot_galex_nuv()
            elif i == 2:
                try:
                    self.plot_cs()
                except:
                    print('no cs halpha')
            elif (i == 3):
                wband = i+1
                self.plot_unwise(band=wband)
            if i > 0:
                plt.gca().set_yticks(())
            #plt.gca().set_xticks(())            
        plt.text(-1.,-.1,self.rootname,transform=plt.gca().transAxes,fontsize=14,horizontalalignment='center')    
        plt.savefig(self.rootname+'-sfr-cutouts.png')
        plt.savefig(self.rootname+'-sfr-cutouts.pdf')        
        
    def plot_legacy_jpg(self):
        # plot jpeg from legacy survey
        t = Image.open(self.legacy_jpegname)        
        plt.imshow(t,origin='upper')
        plt.title(r'$Legacy$')        
        pass
    def plot_legacy(self,band=1):
        # plot image from legacy survey
        # band refers to grz (1,2,3)
        if band == 1:
            display_image(self.legacy_g)
            plt.title(r'$Legacy \ g$')            
        elif band == 2:
            display_image(self.legacy_r)
            plt.title(r'$Legacy \ r$')            
        elif band == 3:
            display_image(self.legacy_z)
            plt.title(r'$Legacy \ z$')
        pass
    def plot_unwise(self,band=1):
        # plot image from legacy survey
        # band refers to W1, W2, W3, W4 (1,2,3,4)
        if band == 1:
            display_image(self.w1,lowrange=False)
            plt.title(r'$unWISE \ W1$')            
        elif band == 2:
            display_image(self.w2,lowrange=False)
            plt.title(r'$unWISE \ W2$')            
        elif band == 3:
            display_image(self.w3,lowrange=False)
            plt.title(r'$unWISE \ W3$')
        elif band == 4:
            display_image(self.w4,lowrange=False)
            plt.title(r'$unWISE \ W4$')
        pass
    def plot_ha(self):
        #v1,v2=scoreatpercentile(self.ha,[vmin,vmax])#.5,99
        #Halpha plus continuum
        #plt.imshow(self.ha,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
        display_image(self.ha)
        plt.title(r'$H\alpha + cont$',fontsize=14)
        
    def plot_r(self):
        #R
        #v1,v2=scoreatpercentile(self.r,[vmin,vmax])#.5,99        
        #plt.imshow(self.r,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')

        display_image(self.r)
        pointing = self.header['object'].replace('ointing','').replace('-','')
        plt.title(r'$R$'+' ('+str(pointing)+')',fontsize=14)
        
    def plot_cs(self):
        #v1,v2=scoreatpercentile(self.cs,[vmin,vmax])
        #plt.imshow(self.cs,origin='lower',cmap='gray_r',vmin=v1,vmax=v2)
        #plt.gca().set_yticks(())
        display_image(self.cs,lowrange=False,percent=99.9,sigclip=False)
        plt.title(r'$H\alpha$',fontsize=14)
    def plot_mask(self):
        #v1,v2=scoreatpercentile(self.cs,[vmin,vmax])
        #plt.imshow(self.cs,origin='lower',cmap='gray_r',vmin=v1,vmax=v2)
        #plt.gca().set_yticks(())
        display_image(self.mask)
        plt.title(r'$H\alpha$',fontsize=14)
        
    def plot_galex_nuv(self):
        if self.nuv_flag:
            display_image(self.nuv_image)
        plt.title(r'$GALEX \ NUV$')


    def adap2020(self,plotsingle=True):
        '''cutouts to show in ADAP 2020 proposal, shows legacy color, W1-W4 and galex'''

        nrow = 2
        ncol = 3

        #Continuum subtracted image
        self.legacy_jpegname = 'legacy/VFID0566-NGC5989-legacy-158.jpg'
        self.nuv_image_name = 'galex/VFID0566-NGC5989-NGC5989-nuv-158.fits'        
        self.nuv_image = fits.getdata(self.nuv_image_name)
        titles = [r'$Legacy \ grz$',r'$GALEX \ NUV$',r'$unWISE \ W1$',r'$unWISE \ W2$',r'$unWISE \ W3$',r'$unWISE \ W4$']        
        if plotsingle:
            figure_size=(9,7)
            plt.figure(figsize=figure_size)
            plt.clf()
            plt.subplots_adjust(left=.05,right=.95,bottom=.05,top=.9,hspace=.2,wspace=0)
        for i in range(nrow*ncol):
            
            plt.subplot(nrow,ncol,i+1)
            if i == 0:
                self.plot_legacy_jpg()
            elif i == 1:
                if self.nuv_flag:
                    self.plot_galex_nuv()
            else: 
                wband = i-1
                self.plot_unwise(band=wband)
            x1,x2 = plt.xlim()
            # zoom in by a factor of 2
            dx = (x2-x1)
            print('image size in unwise = ',dx/2)
            xmean = np.mean([x1,x2])
            xmin = xmean - dx/4
            xmax = xmean + dx/4
            if i > 0:
                plt.axis([xmin,xmax,xmax,xmin])
            else:
                plt.axis([xmin,xmax,xmin,xmax])
            print(x1,x2,xmean,xmin,xmax)
            plt.title(titles[i],fontsize=24)                        
            plt.xticks([])
            plt.yticks([])            
        #plt.text(-1.3,3.8,self.galid,transform=plt.gca().transAxes,fontsize=14,horizontalalignment='center')    
        plt.savefig(self.galid+'-adap2020-cutouts.png')
        plt.savefig(self.galid+'-adap2020-cutouts.pdf')        
        
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='This program will create a plot of halpha image.  will also download images from galex, legacy survey, and unwise.', formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--r', dest = 'r', default = None, help = 'R-band cutout image.  This is all you need to provide if images are named: blah_R.fits, blah_CS.fits, blah_Ha.fits')
    parser.add_argument('--plotall', dest = 'plotall', default = False, action='store_true', help = 'set this to download images and generate all three plots (Ha, all images, sfr indicators)')
    #parser.add_argument('--plotall', dest = 'plotall', default = False, action='store_true', help = 'set this to download images and generate all three plots (Ha, all images, sfr indicators)')        
    args = parser.parse_args()

    if args.r is None:
        print('you must supply the r-band image name!')
        print('try again')
        sys.exit()
    c = cutouts(args.r)
    
    if args.plotall:
        c.runall()
        c.plotcutouts()
        c.plotallcutouts()
        c.plotsfrcutouts()    
        c.plotmask()
    else:
        c.get_halpha_cutouts()
        c.plotcutouts()
