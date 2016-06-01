#!/usr/bin/env python

from astropy.io import fits
import numpy as np
from pylab import *
import sys
import os
mosaic_pixelscale=.423

mypath=os.getcwd()
if mypath.find('Users') > -1:
    print "Running on Rose's mac pro"
    homedir='/Users/rfinn/'
elif mypath.find('home') > -1:
    print "Running on coma"
    homedir='/home/share/'
ratio_r2ha=18.43

def plotintens(data_file,pls='-',pcolor='k',pixelscale=.423,scale=1):

    data1 = np.genfromtxt(data_file)
    sma_pix = data1[:,1]*pixelscale
    intens, intens_err = data1[:,2], data1[:,3]
    plot(sma_pix,intens*scale,ls=pls,color=pcolor)
    errorbar(sma_pix,intens*scale,intens_err,fmt=None,ecolor=pcolor)
    axhline(y=0,ls='--',color='k')


def plotfenclosed(data_file,pls='-',pcolor='k',pixelscale=.423):
    data1 = np.genfromtxt(data_file)
    sma_pix = data1[:,1]*pixelscale
    tflux_e = data1[:,21]
    plot(sma_pix,tflux_e/max(tflux_e)*100.,ls=pls,color=pcolor)
    #errorbar(sma_pix,intens,intens_err,fmt=None,ecolor=pcolor)
    axhline(y=0,ls='--',color='k')


def plotisomag(data_file,pls='-',pcolor='k',pixelscale=.423):
    data1 = np.genfromtxt(data_file)
    sma_pix = data1[:,1]*pixelscale
    mag = data1[:,18]
    magerr_l = data1[:,19]
    magerr_u = data1[:,19]
    magerr=array(zip(magerr_l,magerr_u),'f')
    
    plot(sma_pix,mag,ls=pls,color=pcolor)
    errorbar(sma_pix,mag,yerr=magerr.T,fmt=None,ecolor=pcolor)
    axhline(y=0,ls='--',color='k')


def plotimage(fits_image,vmin=0,vmax=4):
    rfits=fits.open(fits_image)
    im=rfits[0].data.copy()
    rfits.close()
    ax=gca()
    axis('equal')
    logscale=0
    if logscale:
        vmax=log10(np.amax(im)*.9)
        print vmax
        im=log10(im)
        vmin=0.
    else:
        print 'image max = ',np.amax(im)
        vmax=np.amax(im)*.005
        vmin=-1
    imshow((im),interpolation='nearest',origin='upper',cmap='binary',vmin=vmin,vmax=vmax)        
    ax.set_yticklabels(([]))
    ax.set_xticklabels(([]))

def putlabel(s):
    text(.08,.9,s,fontsize=16,transform=gca().transAxes,horizontalalignment='left')
    
def makeplots(rimage,haimage):
    # subplot dimensions
    nx=4
    ny=1

    figure(figsize=(12,3.5))
    subplots_adjust(left=.02,right=.95,wspace=.25,bottom=.2,top=.9)

    subplot(ny,nx,1)

    # plot r-band cutout image
    rfits=rimage
    plotimage(rfits,vmin=-.05,vmax=500)
    putlabel('$R-band $')
    t=rimage.split('_')
    agcnumber=t[0]
    xlabel('$'+str(agcnumber)+'$',fontsize=20)
    subplot(ny,nx,2)
    # plot 24um cutout
    hafits=haimage
    plotimage(hafits,vmin=-.05,vmax=100)
    putlabel(r'$H \alpha $')
    
    subplot(ny,nx,3)
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
    gca().set_yscale('log')
    #gca().set_xscale('log')
    axis([.3,100,.1,1000])
    #axis([.3,40,-10,400.])
    xlabel('$ sma \ (arcsec) $',fontsize=20)
    putlabel('$Intensity $')
    # 21 total flux enclosed by ellipse
    subplot(ny,nx,4)

    plotfenclosed(rdat,pixelscale=mosaic_pixelscale,pcolor='b')
    plotfenclosed(hadat,pixelscale=mosaic_pixelscale,pcolor='r')
    try:
        plotfenclosed(chadat,pixelscale=mosaic_pixelscale,pcolor='m')
    except:
        print 'no convolved Halpha'

    ax=gca()
    axis([.3,100,5.,120.])
    ax.set_yscale('log')
    #ax.set_xscale('log')
    axis([.3,40,5.,120.])

    axhline(y=50,ls=':',color='k')
    axhline(y=70,ls=':',color='k')
    xlabel('$ sma \ (arcsec) $',fontsize=20)
    putlabel('$Flux(<r) $')

    outfile=homedir+'research/LocalClusters/Halpha/EllipseProfiles/'+str(agcnumber)+'-ellipse-profiles.png'
    print outfile
    savefig(outfile)
    close()


#haimage=sys.argv[2]
#rimage=sys.argv[1]
agcnumber=sys.argv[1]
rimage=agcnumber+'_R.fits'
haimage=agcnumber+'_Ha.fits'


makeplots(rimage,haimage)
