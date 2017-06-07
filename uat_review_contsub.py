#!/usr/bin/env python


'''

GOAL:

The goal of this program is to have the user identify galaxies
with Halpha emission and to move the cutouts of these galaxies
to a designated directory for further analysis.


PROCEDURE:

- create a list of all continuum-subtracted images in the current directory.
- Display R, Halpha+cont and continuum-subtracted images for each galaxy, one at a time.
- ask user to determine if galaxy has Halpha emission.
- move cutouts of galaxies with Halpha emission to a designated directory for further analysis

EXAMPLE:

if necessary, change file to executable:

chmod +x ~/github/HalphaImaging/uat_review_contsub.py

to run, type the following in the command line:

~/github/HalphaImaging/uat_review_contsub.py --cluster 'A1367'

To get cutouts from a particular pointing:

~/github/HalphaImaging/uat_review_contsub.py --cluster 'A1367-h02'

or if you want to put cutouts of Halpha emitters in a central directory:

~/github/HalphaImaging/uat_review_contsub.py --cluster 'A1367' --dir '/home/share/research/HalphaImaging/HalphaEmitters'

INPUT/OUPUT:

REQUIRED MODULES:

EXTRA NOTES:

WRITTEN BY:

Rose Finn

'''
from astropy.io import fits
import argparse
from matplotlib import pyplot as plt
from scipy.stats import scoreatpercentile
import glob
import os
import errno
import sys

parser = argparse.ArgumentParser(description ='Review continuum-subtracted images to determine if galaxy has Halpha emission.')
parser.add_argument('--cluster', dest = 'cluster', default = None, help = 'cluster and prefix of image names (e.g. A1367 or A1367-h02)')
parser.add_argument('--dir', dest = 'output_dir', default = 'HalphaEmitters', help = 'Directory to move Halpha emitters to.  Default is a subdirectory of the working directory called HalphaEmitters.')
args = parser.parse_args()

try:
    os.makedirs(args.output_dir)
except OSError as exception:
    if exception.errno != errno.EEXIST:
        raise
    else:
        print 'Just to let you know, the output directory already exists.'
        print 'If this is not what you want, remove the output directory and start over.\n'

    
figure_size=(10,4)

# gather the list of continuum-subtracted images
infiles = glob.glob(args.cluster+'*-CS.fits')

for csimage in infiles:
    prefix = csimage.split('-CS.fits')[0]
    print 'prefix = ',prefix
    haimage = prefix+'*-Ha.fits'
    rimage = prefix+'*-R.fits'

    r,r_header = fits.getdata(rimage,header=True)
    ha,ha_header = fits.getdata(haimage,header=True)
    cs,cs_header = fits.getdata(csimage,header=True)
    v1,v2=scoreatpercentile(ha,[5.,95.])
    plt.figure(1,figsize=figure_size)
    plt.clf()
    plt.subplots_adjust(hspace=0,wspace=0)
    #Halpha plus continuum
    plt.subplot(1,3,1)
    plt.imshow(ha,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
    plt.title('Halpha + cont')
    plt.xlabel(prefix, fontsize = 14)
    #R
    plt.subplot(1,3,2)
    plt.imshow(r,cmap='gray_r',vmin=v1/.0445,vmax=v2/.0445,origin='lower')
    plt.title('R')
    plt.gca().set_yticks(())
    #Continuum subtracted image
    plt.subplot(1,3,3)
    plt.imshow(cs,cmap='gray_r',vmin=v1,vmax=v2,origin='lower')
    plt.title('contsub')
    plt.gca().set_yticks(())
    plt.show(block=False)
    need_input = True
    while need_input:
        t=raw_input('Does continuum-subtracted image show Halpha emission? (y or n, q to quit) \n')
        if t.find('n') > -1: # no Halpha emission
            need_input = False
            plt.savefig(prefix+'-noHa.png')
        elif t.find('y') > -1:
            # move cutouts to specified directory
            need_input = False
            os.rename(rimage,args.output_dir+'/'+rimage)
            os.rename(haimage,args.output_dir+'/'+haimage)
            os.rename(csimage,args.output_dir+'/'+csimage)
            plt.savefig(args.output_dir+'/'+prefix+'-Ha-emission.png')
        elif t.find('q') > -1:
            need_input = False
            sys.exit()
        else:
            print "sorry, I didn't understand what you said."
            print "let's try again...\n"
            
