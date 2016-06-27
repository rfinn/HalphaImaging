#!/usr/bin/env python

'''

GOAL:

PROCEDURE:

EXAMPLE:

from within ipython -pylab

%run ~/github/HalphaImaging/uat_cutouts_with_ellipse_profiles.py --r A1367-140231-R.fits --ha A1367-140231-Ha.fits

REQUIRED MODULES:

NOTES:

WRITTEN BY:
Rose Finn
updated 2016-06-23

'''


from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import argparse
from scipy.stats import scoreatpercentile
mosaic_pixelscale=.425


parser = argparse.ArgumentParser(description ='Run ellipse on R and Halpha images')
#parser.add_argument('--r', dest = 'r', default = None, help = 'R-band image')
#parser.add_argument('--ha', dest = 'ha', default = None, help = 'R-band image')
parser.add_argument('--cluster', dest = 'cluster', default = None, help = 'cluster and prefix of image names (e.g. A1367)')
parser.add_argument('--id', dest = 'id', default = None, help = 'NSA ID number')

#parser.add_argument('--scale', dest = 'scale', default = 0.06, help = 'R-band image')
args = parser.parse_args()


ratio_r2ha=18.43

def plotintens(data_file,pls='-',pcolor='k',pixelscale=.423,scale=1):

    data1 = np.genfromtxt(data_file)
    sma_pix = data1[:,1]*pixelscale
    intens, intens_err = data1[:,2], data1[:,3]
    plt.plot(sma_pix,intens*scale,ls=pls,color=pcolor)
    plt.errorbar(sma_pix,intens*scale,intens_err,fmt=None,ecolor=pcolor)
    plt.axhline(y=0,ls='--',color='k')


def plotfenclosed(data_file,pls='-',pcolor='k',pixelscale=0.425):
    data1 = np.genfromtxt(data_file)
    sma_pix = data1[:,1]*pixelscale
    tflux_e = data1[:,21]
    plt.plot(sma_pix,tflux_e/max(tflux_e)*100.,ls=pls,color=pcolor)
    #errorbar(sma_pix,intens,intens_err,fmt=None,ecolor=pcolor)
    plt.axhline(y=0,ls='--',color='k')


def plotisomag(data_file,pls='-',pcolor='k',pixelscale=.423):
    data1 = np.genfromtxt(data_file)
    sma_pix = data1[:,1]*pixelscale
    mag = data1[:,18]
    magerr_l = data1[:,19]
    magerr_u = data1[:,19]
    magerr=array(zip(magerr_l,magerr_u),'f')
    
    plt.plot(sma_pix,mag,ls=pls,color=pcolor)
    plt.errorbar(sma_pix,mag,yerr=magerr.T,fmt=None,ecolor=pcolor)
    plt.axhline(y=0,ls='--',color='k')


def plotimage(fits_image,vmin=0,vmax=4):
    im=fits.getdata(fits_image)
    ax=plt.gca()
    plt.axis('equal')
    logscale=0
    vmin,vmax=scoreatpercentile(im,[5.,95.])
    plt.imshow((im),interpolation='nearest',origin='upper',cmap='binary',vmin=vmin,vmax=vmax)        
    ax.set_yticklabels(([]))
    ax.set_xticklabels(([]))

def putlabel(s):
    plt.text(.08,.9,s,fontsize=16,transform=plt.gca().transAxes,horizontalalignment='left')
    
def makeplots(rimage,haimage):
    # subplot dimensions
    nx=4
    ny=1

    plt.figure(figsize=(12,3.5))
    plt.subplots_adjust(left=.02,right=.95,wspace=.25,bottom=.2,top=.9)

    plt.subplot(ny,nx,1)

    # plot r-band cutout image
    rfits=rimage
    plotimage(rfits,vmin=-.05,vmax=500)
    putlabel('$R-band $')
    t=rimage.split('-')
    agcnumber=t[1]
    plt.xlabel('$'+str(agcnumber)+'$',fontsize=20)
    plt.subplot(ny,nx,2)
    # plot 24um cutout
    hafits=haimage
    plotimage(hafits,vmin=-.05,vmax=100)
    putlabel(r'$H \alpha $')
    
    plt.subplot(ny,nx,3)
    # plot r and 24 profiles
    t=rimage.split('fits')
    rdat=t[0]+'dat'

    t=haimage.split('fits')
    hadat=t[0]+'dat'

    chadat='c'+hadat

    plotintens(rdat,pixelscale=mosaic_pixelscale,pcolor='b')
    plotintens(hadat,pixelscale=mosaic_pixelscale,pcolor='r',scale=ratio_r2ha)

    try:
        plotintens(chadat,pixelscale=mosaic_pixelscale,pcolor='m',scale=ratio_r2ha)
    except:
        print 'no convolved Halpha'
    plt.gca().set_yscale('log')
    #gca().set_xscale('log')
    #plt.axis([.3,100,.1,1000])
    #axis([.3,40,-10,400.])
    plt.xlabel('$ sma \ (arcsec) $',fontsize=20)
    putlabel('$Intensity $')
    # 21 total flux enclosed by ellipse
    plt.subplot(ny,nx,4)

    plotfenclosed(rdat,pixelscale=mosaic_pixelscale,pcolor='b')
    plotfenclosed(hadat,pixelscale=mosaic_pixelscale,pcolor='r')
    try:
        plotfenclosed(chadat,pixelscale=mosaic_pixelscale,pcolor='m')
    except:
        print 'no convolved Halpha'

    ax=plt.gca()
    plt.axis([.3,10,5.,120.])
    ax.set_yscale('log')
    #ax.set_xscale('log')
    plt.axis([.3,40,5.,120.])

    plt.axhline(y=50,ls=':',color='k')
    plt.axhline(y=70,ls=':',color='k')
    plt.xlabel('$ sma \ (arcsec) $',fontsize=20)
    putlabel('$Flux(<r) $')

    outfile=prefix+'-ellipse-profiles.png'
    print 'saving result as ',outfile
    plt.savefig(outfile)
    #plt.close()


rimage=args.cluster+'-'+args.id+'-R.fits'
haimage=args.cluster+'-'+args.id+'-CS.fits'
#t = rimage.split('-')
prefix = args.cluster#t[0]+'-'+t[1] # should be CLUSTER-NSAID
#haimage=args.ha


makeplots(rimage,haimage)
