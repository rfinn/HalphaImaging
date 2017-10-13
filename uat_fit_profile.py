#!/usr/bin/env python

'''
GOAL:
   - To measure the effective radius on a galaxy based on its intensity vs. radius

PROCEDURE:
   - Read in photometry files from uat_measure_ellip_phot.py
   - Read in radii from SExtractor catalog
   - Fit a linear function to log10(intensity) vs radius in arcsec

BACKGROUND:

EXAMPLE:

from within ipython

%run ~/github/HalphaImaging/uat_fit_profile.py --cluster pointing-1 --id 118647 --compareNSA

WRITTEN BY:
   Rose A. Finn
   July 06, 2017



'''
import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
import argparse
import scipy
import scipy.optimize

def simple_fit(plotsingle=False):
    # function to test fitting in log space
    data1 = np.genfromtxt('pointing-1-118647-R_phot.dat')
    linear_fit(data1[:,0],np.log(data1[:,2]),plotsingle=plotsingle)
    
def linear_fit(x,y,plotsingle=False):
    # testing linear fit to log(intensity) vs radius
    if plotsingle:
        plt.figure()
    cc=np.polyfit(x,y,1)
    xl=np.linspace(0,np.max(x),100)
    yl=np.polyval(cc,xl)
    plt.plot(x,y)
    plt.plot(xl,yl,'r-')
    return cc,xl,yl

def sersic(r,I0,n,k):
    return I0*np.exp(-1*k*r**(1./n))

def sersic_neq1(r,I0,k):
    return I0*np.exp(-k*r)

def read_phot_file(pfile):
    # read in file from uat_measure_ellip_phot
    # return radius in pixels and intensity
    dat = np.genfromtxt(pfile)
    return dat[:,0],dat[:,2]

parser = argparse.ArgumentParser(description ='Fit a Sersic function to R and Halpha radial profiles')
parser.add_argument('--cluster', dest = 'cluster', default = None, help = 'cluster and prefix of image names (e.g. A1367)')
parser.add_argument('--id', dest = 'id', default = None, help = 'NSA ID number')
parser.add_argument('--scale', dest = 'scale', default = 0.0445, help = 'Filter ratio of R-band to Halpha.  Default is 0.0445.')
parser.add_argument('--pixelscale', dest = 'pixelscale', default = 0.43, help = 'Pixel scale.  Default is set to HDI pixel scale of 0.43 arcsec/pixel.')
parser.add_argument('--rmax', dest = 'rmax', default = 4., help = 'maximum radius to fit profile out to.  This is a multiple of what sextractor measures as R90.  Default is 4xR90.')
parser.add_argument('--compareNSA', dest = 'compareNSA', default = False, action='store_true', help = 'Compare results to the NSA fit')
parser.add_argument('--NSApath', dest = 'NSApath', default = '/Users/rfinn/research/NSA/',  help = 'Path to NSA file.  Default is ~/research/NSA/.  Set this if using compareNSA.')
args = parser.parse_args()


if __name__ == '__main__':
        
    # this is the main part of the program

    # names of R-band and continuum-subtracted Halpha photometry files
    filenameR = args.cluster+'-'+args.id+'-R_phot.dat'
    filenameHa = args.cluster+'-'+args.id+'-CS_phot.dat'

    # read in R-band SE cat
    secatR =  args.cluster+'-'+args.id+'-R.cat'
    cat = fits.getdata(secatR,2)
    # read in Halpha SE cat
    secat_ha =  args.cluster+'-'+args.id+'-CS.cat'
    cat_ha = fits.getdata(secat_ha,2)

    # get dimensions of image
    image =  args.cluster+'-'+args.id+'-R.fits'
    imdat = fits.getdata(image)
    xdim,ydim = imdat.shape
    
    # locate central object in sextractor catalog
    distance = np.sqrt((cat.X_IMAGE - xdim/2.)**2 + (cat.Y_IMAGE - ydim/2.)**2)
    objectID = distance == min(distance)

    # find R90 - radius enclosing 90% of flux

    R90 = cat.FLUX_RADIUS[objectID][0][0]*args.pixelscale
    R50 = cat.FLUX_RADIUS[objectID][0][1]*args.pixelscale
    R30 = cat.FLUX_RADIUS[objectID][0][2]*args.pixelscale

    R90_ha = cat.FLUX_RADIUS[objectID][0][0]*args.pixelscale
    R50_ha = cat.FLUX_RADIUS[objectID][0][1]*args.pixelscale
    R30_ha = cat.FLUX_RADIUS[objectID][0][2]*args.pixelscale

    
    # set maximum for fitting radius
    rmax = args.rmax*R90  # selected this to give a Re that was relatively close to NSA value

    # read in Rband file
    radius,intensity = read_phot_file(filenameR)
    intensity = intensity/max(intensity)
    
    # read in CS Halpha file
    radius_ha, intensity_ha = read_phot_file(filenameHa)
    intensity_ha = intensity_ha/max(intensity_ha)
    
    # convert radius to arcsec
    radius = radius*args.pixelscale
    radius_ha = radius_ha*args.pixelscale

        
    # limit radius for fitting to < 2*R90, based on Rband image
    flag = (radius > 0.25*rmax) & (radius < 3.5*rmax)
    flag = radius < rmax

    popt, pcov = scipy.optimize.curve_fit(sersic_neq1,radius[flag],intensity[flag])
    popt_ha, pcov_ha = scipy.optimize.curve_fit(sersic_neq1,radius_ha[flag],intensity_ha[flag])

    # assume exponential and do linear fit to log(intensity) vs. radius
    cc = np.polyfit(radius[flag],np.log(intensity[flag]),1)
    flag = flag & (intensity_ha > 0)
    cc_ha = np.polyfit(radius_ha[flag],np.log(intensity_ha[flag]),1)

    xl=np.linspace(0,np.max(radius),100)
    yl = np.polyval(cc,xl)
    yl_ha = np.polyval(cc_ha,xl)
    
    
    plt.figure(figsize=(6,8))
    plt.subplots_adjust(top=.925)
    
    for i in range(3):
        plt.subplot(3,1,i+1)
        if i == 0:
            plt.plot(radius, intensity, 'b.', label='R',markersize=6)
            plt.plot(radius,sersic_neq1(radius, *popt),'r-',label='exp fit')
            #c,x,y = linear_fit(radius,np.log(intensi
            plt.plot(xl,np.exp(yl),'r:',label='linear fit to log(I)')
            plt.axvline(x=R50,ls='--',label='SE R50')
            plt.axvline(x=-1./cc[0],ls=':',label='log fit R50')
            plt.axvline(x=1./popt[1],ls='-.',label='exp fit R50')

        elif i == 1:
            plt.plot(radius_ha, intensity_ha, 'c.', label='Ha',markersize=6)
            plt.plot(radius_ha,sersic_neq1(radius_ha, *popt_ha),'r-',label='exp fit')
            plt.plot(xl,np.exp(yl_ha),'r:',label='linear fit to log(I)')
            plt.axvline(x=R50_ha,ls='--',label='SE R50')
            plt.axvline(x=-1./cc_ha[0],ls=':',label='log fit R50')
            plt.axvline(x=1./popt_ha[1],ls='-.',label='exp fit R50')
        elif i == 2:
            # plot R and Halpha on the same axis
            plt.plot(radius, intensity, 'b.', label='R',markersize=6)
            plt.plot(radius_ha, intensity_ha, 'c.', label='Ha',markersize=6)
        plt.axvline(x=rmax,color='k',ls='-',label='max r for fit')
        plt.legend(loc = 'upper right',prop={'size':7})

        plt.gca().set_yscale('log')
        if i == 0:
            #plt.ylabel('$ R-band \ Intensity \ (ADU/pixel^2)$')
            plt.title('NSA ID:'+args.id)
        elif i == 1:
            plt.ylabel('$Normalized \ Intensity \ (ADU/pixel^2)$')
        elif i == 2:
            plt.xlabel('Radius (arcsec)')
        plt.ylim(1.e-4,2)
        plt.xlim(-2,np.max(radius_ha))
    plt.savefig(args.cluster+'-'+args.id+'-radial-profile.png')
    print 'R-BAND FIT RESULTS:'
    print 'I0 = %.2f'%(popt[0])
    #print 'n =  %.2f'%(popt[1])
    print 'Re (arcsec)= %.2f'%(1./popt[1])

    print'\nLINEAR FIT TO LOG(I)'
    print 'I0 = %.2f'%(np.exp(cc[1]))
    #print 'n =  %.2f'%(popt[1])
    print 'Re (arcsec)= %.2f'%(-1./cc[0])

    print '\nSExtractor Estimate for R50'
    print 'Re (arcsec)= %.2f'%(R50)
    
    if args.compareNSA:
        # read in NSA file
        nsafile = args.NSApath+'nsa_v0_1_2.fits'
        nsa = fits.getdata(nsafile)
        nsa_index = np.where(nsa.NSAID == int(args.id))
        print '\nNSA fit results'

        print 'n = %.2f'%(nsa.SERSIC_N[nsa_index][0])
        print 'theta 50 (arcsec)= %.2f'%(nsa.SERSIC_TH50[nsa_index][0])
    print '\nH-ALPHA FIT RESULTS:'
    print 'I0 = %.2f'%(popt_ha[0])
    #print 'n =  %.2f'%(popt_ha[1])
    print 'Re (arcsec)= %.2f'%(1./popt_ha[1])


    print'\nLINEAR FIT TO LOG(I)'
    print 'I0 = %.2f'%(np.exp(cc_ha[1]))
    #print 'n =  %.2f'%(popt[1])
    print 'Re (arcsec)= %.2f'%(-1./cc_ha[0])

    print '\nSExtractor Estimate for R50'
    print 'Re (arcsec)= %.2f'%(R50_ha)
