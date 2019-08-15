#!/usr/bin/env python

'''
PURPOSE:

The goal of the program is to create a mask for a galaxy image to mask
out other objects within the cutout area.

USAGE:

from within ipython type:

   %run ~/github/HalphaImaging/uat_mask.py --image 'A1367-140231-R.fits'

you just need to run this on R-band images.

Interacting with the display is finicky.  This works fine when running
within ipython - not so much when running from the command line.  When running
from the command line, I am not able to interact with the figure.  This may
have something to do with setting block=False in show().


PROCEDURE:


REQUIRED MODULES:
   os
   astropy
   numpy
   argsparse
   matplotlib
   scipy

NOTES:
- rewrote using a class

'''

import os
import sys
from astropy.io import fits
import numpy as np
import argparse
import pyds9
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile

defaultcat='default.sex.HDI.mask'


class mask_image():
    def __init__(self, image=None, haimage=None, sepath='~/github/HalphaImaging/astromatic/', nods9=False,
                 param='default.sex.HDI.mask', threshold=0.05,snr=2,cmap='gist_heat_r'):

        self.image_name = image
        if haimage != None:
            self.haimage_name = haimage
        self.sepath = sepath
        self.nods9 = nods9
        self.param = param
        self.threshold = threshold
        self.snr = snr        
        self.cmap = cmap
        self.xcursor_old = -99
        self.xcursor = -99

        # create name for output mask file
        t = self.image_name.split('.fit')
        self.mask_image=t[0]+'-mask.fits'
        self.mask_inv_image=t[0]+'-inv-mask.fits'
        print('saving image as: ',self.mask_image)
        
        # read in image and define center coords
        self.image, self.imheader = fits.getdata(self.image_name,header = True)
        self.ymax,self.xmax = self.image.shape
        self.xc = self.xmax/2.
        self.yc = self.ymax/2.

        self.v1,self.v2=scoreatpercentile(self.image,[5.,99.5])
        self.adjust_mask = True
        self.figure_size = (10,5)
        self.mask_size = 20. # side of square to mask out when user clicks on a pixel

        # set up array to store the user-created object masks
        self.usr_mask = np.zeros_like(self.image)

        # set off center flag as false by default
        self.off_center_flag = False

        # keep track of extra objects that the user deletes from mask

        self.deleted_objects = []
        self.link_files()
        if not(self.nods9):
             # open ds9
            try:
                d.set('frame delete all')
            except NameError:
                d=pyds9.DS9()
                d.set('frame delete all')


    def link_files(self):
        sextractor_files=['default.sex.HDI.mask','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('ln -s '+self.sepath+'/'+file+' .')
    def clean_links(self):
        # clean up
        #sextractor_files=['default.sex.sdss','default.param','default.conv','default.nnw']
        sextractor_files=['default.sex.HDI.mask','default.param','default.conv','default.nnw']
        for file in sextractor_files:
            os.system('unlink '+file)

    def read_se_cat(self):
        sexout=fits.getdata('test.cat')
        self.xsex=sexout['XWIN_IMAGE']
        self.ysex=sexout['YWIN_IMAGE']
        self.fwhm = sexout['FWHM_IMAGE']
        dist=np.sqrt((self.yc-self.ysex)**2+(self.xc-self.xsex)**2)
        #   find object ID
        objIndex=np.where(dist == min(dist))
        objNumber=sexout['NUMBER'][objIndex]
        return objNumber[0] # not sure why above line returns a list

    def runse(self,galaxy_id = None):
        print('using a deblending threshold = ',self.threshold)
        os.system('sex %s -c %s -CATALOG_NAME test.cat -CATALOG_TYPE FITS_1.0 -DEBLEND_MINCONT %f -DETECT_THRESH %f -ANALYSIS_THRESH %f'%(self.image_name,self.param,float(self.threshold),float(self.snr),float(self.snr)))
        self.maskdat = fits.getdata('segmentation.fits')
        # grow masked areas
        bool_array = np.array(self.maskdat.shape,'bool')
        #for i in range(len(self.xsex)):
            
            
        if self.off_center_flag:
            print('setting center object to objid ',self.galaxy_id)
            self.center_object = self.galaxy_id
        else:
            self.center_object = self.read_se_cat()
        self.maskdat[self.maskdat == self.center_object] = 0
        # add back the square masked areas that the user added
        self.maskdat = self.maskdat + self.usr_mask
        # remove objects that have already been deleted by user
        if len(self.deleted_objects) > 0:
            for objID in self.deleted_objects:
                self.maskdat[self.maskdat == objID] = 0.
        # write out mask image
        fits.writeto(self.mask_image,self.maskdat,header = self.imheader,overwrite=True)
        invmask = self.maskdat > 0.
        invmask = np.array(~invmask,'i')
        fits.writeto(self.mask_inv_image,invmask,header = self.imheader,overwrite=True)
    def show_mask(self):
        if args.nods9:
            plt.close('all')
            self.fig = plt.figure(1,figsize=self.figure_size)
            plt.clf()
            plt.subplots_adjust(hspace=0,wspace=0)
            plt.subplot(1,2,1)
            plt.imshow(self.image,cmap='gray_r',vmin=self.v1,vmax=self.v2,origin='lower')
            plt.title('image')
            plt.subplot(1,2,2)
            #plt.imshow(maskdat,cmap='gray_r',origin='lower')
            plt.imshow(self.maskdat,cmap=args.cmap,origin='lower')
            plt.title('mask')
            plt.gca().set_yticks(())
            #plt.draw()
            plt.show(block=False)
        else:
            self.ds9_open()
            s='file new '+self.image_name
            self.d.set(s)
            self.ds9_adjust()

            if args.haimage != None:
                s='file new '+self.haimage_name
                self.d.set(s)
                self.ds9_adjust()
            s='file new '+self.mask_image
            self.d.set(s)
            self.ds9_adjust()
 
        
    def ds9_adjust(self):
        self.d.set('scale log')
        self.d.set('zoom to fit')
        self.d.set('scale mode 99.5')      

    def ds9_onclick(self):
        s=self.d.get('iexam')
        print(s)
        a,b = s.split()
        print(a)
        print(b)
        self.xcursor = int(float(a))
        self.ycursor = int(float(b))
    def onclick(self,event):
        print('xdata=%f, ydata=%f' %(event.xdata, event.ydata))
        self.xcursor = event.xdata
        self.ycursor =  event.ydata

    def get_usr_mask(self):
        print('click on the location to add object mask')
        if args.nods9:
            a = self.fig.canvas.mpl_connect('button_press_event', self.onclick)
    
            while self.xcursor == self.xcursor_old: # stay in while loop until mouse is clicked
                plt.pause(1)
        else:
            self.ds9_onclick()
        # mask out a rectangle around click
        # size is given by mask_size
        xmin = int(self.xcursor) - int(0.5*self.mask_size)
        ymin = int(self.ycursor) - int(0.5*self.mask_size)
        xmax = int(self.xcursor) + int(0.5*self.mask_size)
        ymax = int(self.ycursor) + int(0.5*self.mask_size)
        # make sure mask dimensions are not outside of the image
        if xmin < 0:
            xmin = 0
        if ymin < 0:
            ymin = 0
        if xmax > self.xmax:
            xmax = self.xmax
        if ymax > self.ymax:
            ymax = self.ymax
        print('xcursor, ycursor = ',self.xcursor, self.ycursor)
        mask_value = np.max(self.maskdat) + 1
        self.usr_mask[ymin:ymin+int(self.mask_size),xmin:xmin+int(self.mask_size)] = mask_value*np.ones([self.mask_size,self.mask_size])
        self.maskdat = self.maskdat + self.usr_mask
        fits.writeto(self.mask_image, self.maskdat, header = self.imheader, overwrite=True)
        print('added mask object '+str(mask_value))
        self.xcursor_old = self.xcursor 

    def print_menu(self):
        t=input('enter:\n \t pixel value to remove object in mask;\n \t o if target is off center (and program is removing the wrong object);\n \t a to mask additional pixels;\n \t r to change the size of the mask box;\n \t t to adjust SE threshold (0=lots, 1=no deblend ); \n \t s to adjust SE SNR; \n \t w to write output and quit; \n \t q to quit without saving\n')
        try:
            objID = int(t)
            self.maskdat[self.maskdat == objID] = 0.
            self.deleted_objects.append(objID)
        except ValueError:
            adjust_scale = False
            if t.find('t') > -1:
                t = raw_input('enter new threshold')
                args.threshold = float(t)
                self.runse()
            if t.find('s') > -1:
                t = raw_input('enter new SNR')
                self.snr = float(t)
                self.runse()
            
            if t.find('a') > -1:

                self.get_usr_mask()
        if t.find('r') > -1:
            print('current box size = '+str(self.mask_size))
            t = raw_input('enter new size for square area to be masked (in pixels)\n')
            try:
                self.mask_size = float(t)
            except:
                print('error reading input')
            
        if t.find('o') > -1:
            t = raw_input('enter object number for target galaxy\n')
            self.off_center_flag = True
            self.galaxy_id = int(t)
            self.runse()
        if t.find('q') > -1:
            sys.exit()
        if t.find('w') > -1:
            newfile = fits.PrimaryHDU()
            newfile.data = self.maskdat
            newfile.header = self.imheader
            fits.writeto(self.mask_image, newfile.data, header = newfile.header, overwrite=True)
            self.adjust_mask = False
    def edit_mask(self):
        self.runse()
        while self.adjust_mask:    
            self.show_mask()
            self.print_menu()
            fits.writeto(self.mask_image,self.maskdat,header = self.imheader,overwrite=True)
            
    def ds9_open(self):
        try:
            self.d.set('frame delete all')
        except AttributeError:
            self.d=pyds9.DS9()
            self.d.set('frame delete all')
                
   

        
# run sextractor on input image
# return segmentation image with central object removed

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description ='Create a mask for extraneous objects in field')
    parser.add_argument('--R', dest = 'image', default = None, help = 'R-band image to mask')
    parser.add_argument('--Ha', dest = 'haimage', default = None, help = 'Halpha image to mask to use as a comparison when identifying stars')
    parser.add_argument('--path',dest = 'path', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of default config files.  Default is ~/github/HalphaImaging/astromatic')
    parser.add_argument('--nods9',dest = 'nods9', default = False, action="store_true", help = 'Set this if you DO NOT want to use ds9 to display mask')
    parser.add_argument('--param',dest = 'param', default ='default.sex.HDI.mask', help = 'sextractor parameter file.  Default is default.sex.HDI.mask')
    parser.add_argument('--threshold', dest = 'threshold', default = .005, help = "sextractor DEBLEND_MINCONT: 0=lots of deblending; 1=none (default = .005)",action="store")
    parser.add_argument('--snr', dest = 'snr', default = 2, help = "snr to use for sextractor detections (default is 2).",action="store")
    parser.add_argument('--cmap', dest = 'cmap', default = 'gist_heat_r', help = "color map to use when displaying image mask.  default is gist_heat_r.") 
    args = parser.parse_args()






    
    m = mask_image(args.image, haimage=args.haimage, sepath=args.path, nods9=args.nods9,
                 param=args.param, threshold=args.threshold,snr=args.snr,cmap=args.cmap)
    m.edit_mask()
    m.clean_links()


    




    
