#!/usr/bin/env python

'''

'''
from astropy.io import ascii
from astropy.io import fits
import numpy as np
from matplotlib import pyplot as plt
import getzp
import argparse

parser = argparse.ArgumentParser(description ='Run sextractor, get Pan-STARRS catalog, and then computer photometric ZP\n \n from within ipython: \n %run ~/github/Virgo/programs/getzp.py --image pointing031-r.coadd.fits --instrument i \n\n then:\n x,y = fitzp() \n \n The y intercept is -1*ZP. \n \n x and y data are returned in case you want to make additional plots.', formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--fitzp', dest = 'fitzp', default = False, action= 'store_true', help = 'set this to refit the ZP for each standard star image')

parser.add_argument('--useri',dest = 'useri', default = False, action = 'store_true', help = 'set this to use r->R transformation as a function of r-i rather than the g-r relation.  g-r is the default.')

parser.add_argument('--mag', dest = 'mag', default = 0,help = "select SE magnitude to use when solving for ZP.  0=MAG_APER,1=MAG_BEST,2=MAG_PETRO.  Default is MAG_APER ")
parser.add_argument('--naper', dest = 'naper', default = 5,help = "select fixed aperture magnitude.  0=10pix,1=12pix,2=15pix,3=20pix,4=25pix,5=30pix.  Default is 5 (30 pixel diameter)")
args = parser.parse_args()
args.mag = int(args.mag)
args.naper = int(args.naper)
# images containing standard stars

standard_images = ['nshdztrc7131t0019o00.fits', 'nshdztrc7131t0020o00.fits','nshdztrc7131t0073o00.fits','nshdztrc7131t0074o00.fits']
secats = ['nshdztrc7131t0019o00.cat','nshdztrc7131t0020o00.cat','nshdztrc7131t0073o00.cat','nshdztrc7131t0074o00.cat']


# default values to use in case ZPs are not measured
myzp = np.array([23.476,23.462,23.462,23.472],'f')
myzperr = 0.1*np.ones(len(standard_images),'f')

# measure ZP
if args.fitzp:
    for i,im in enumerate(standard_images):
        zp = getzp.getzp(im, instrument='h', filter='R', nsigma = 2.5, useri = args.useri, mag = args.mag, naper = args.naper)
        zp.getzp()
        myzp[i] = -1.*zp.zp
        myzperr[i] = zp.zperr
zpdict = dict((a,b) for a,b in zip(standard_images,myzp))
zperrdict = dict((a,b) for a,b in zip(standard_images,myzperr))

              

# read in positions and mags of Landolt standards
stand = ascii.read('mystandards.csv', data_start=1, delimiter=',')
sfile = stand['col1']
landolt_rmag = stand['col5']
ximage = np.array(stand['col6'],'f')
yimage = np.array(stand['col7'],'f')
zp = np.zeros(len(sfile),'f')
zperr = np.zeros(len(sfile),'f')
for i,f in enumerate(sfile):
    zp[i] = zpdict[f]
    zperr[i] = zperrdict[f]


# read in SE cats

# read in one catalogs

cat = fits.getdata(secats[0],2)
    
# create array with nrow = nstandards, and format like SE cat
newarray = np.zeros(len(stand),dtype = cat.dtype)

# loop through SE catalogs


for c in secats:
    cat = fits.getdata(c,2)
    # for appropriate standards
    prefix = c.split('.')[0]
    for i,s in enumerate(sfile):
        if s.find(prefix) > -1:
            # find location of landolt standards in SE cats
            distance = np.sqrt((ximage[i] - cat['X_IMAGE'])**2 + (yimage[i] - cat['Y_IMAGE'])**2)
            imatch = np.where(distance == min(distance))
            # populate new array with matched lines in SE cat
            newarray[i] = cat[imatch]
                              


# plot measured mag vs known mag
if args.mag == 0:

    testmag = newarray['MAG_APER'][:,self.naper]
    testmagerr = newarray['MAG_APER'][:,self.naper]
elif args.mag == 1:
    testmag = newarray['MAG_BEST']
    testmagerr = newarray['MAGERR_BEST']
elif args.mag == 2:
    testmag = newarray['MAG_PETRO']
    testmagerr = newarray['MAGERR_PETRO']

plt.figure()
residual = zp + testmag - landolt_rmag
plt.errorbar(landolt_rmag[0:5], residual[0:5],yerr=testmagerr[0:5],fmt='bo',label='PG0918+029')
plt.errorbar(landolt_rmag[5:13], residual[5:13],yerr=testmagerr[5:13],fmt='go',label='RU_149A')
plt.errorbar(landolt_rmag[13:16], residual[13:16],yerr=testmagerr[13:16],fmt='co',label='PG1528')
plt.errorbar(landolt_rmag[16:-1], residual[16:-1],yerr=testmagerr[16:-1],fmt='mo',label='PG1633')
plt.legend()
plt.xlabel('Landolt R mag',fontsize=14)
plt.ylabel('Measured mag',fontsize=14)
plt.axhline(y=0)
plt.figure()
VR = stand['col4']
plt.plot(VR[0:5], residual[0:5],'bo',label='PG0918+029')
plt.plot(VR[5:13], residual[5:13],'go',label='RU_149A')
plt.plot(VR[13:16], residual[13:16],'co',label='PG1528')
plt.plot(VR[16:-1], residual[16:-1],'mo',label='PG1633')
plt.legend()
plt.xlabel('Landolt V-R',fontsize=14)
plt.ylabel('Measured mag - Landolt Mag',fontsize=14)
plt.axhline(y=0)


