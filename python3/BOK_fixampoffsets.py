#!/usr/bin/env python

"""
GOAL:

This is a major rewrite of BOK_pipeline_fixampoffsets.py

to make use of .head file - need to run SE in mef mode

USEAGE: 

python ~/github/HalphaImaging/python3/BOK_fixampoffsets.py imagename 



OUTPUT:
* for each input image, an output image named "m"+imagename will be created.
* the output image is the sky-subtracted image

"""


import sys
import os
import shutil
import numpy as np
from astropy.io import fits
from astropy.table import Table
import astropy.units as u
import astropy.coordinates as coord
from astropy.coordinates import SkyCoord
from astropy.stats import sigma_clip
from astropy.stats import mad_std
from matplotlib import pyplot as plt

homedir = os.getenv("HOME")
sys.path.append(homedir+"/github/HalphaImaging/python3/")
from getzp import panstarrs_query

from scipy.optimize import curve_fit
try:
    from scipy.stats import median_abs_deviation as MAD2
except:
    from scipy.stats import median_absolute_deviation as MAD2

zpfunc = lambda x, zp: x + zp
mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']

def plot_se_pan_positions(shdu,pan,plotprefix="test"):
    """plot positions of SE and panstarrs catalogs  """
    plt.figure(figsize=(10,10))
    plt.plot(pan['RAJ2000'],pan['DEJ2000'],'ko',markersize=10,mfc='None',label='PanSTARRS',alpha=.5)
    icolor = 0
    for i in [2,4,6,8]:
        plt.plot(shdu[i].data['ALPHA_J2000'],shdu[i].data['DELTA_J2000'],'r.',color=mycolors[icolor],label='SE-CCD'+str(i//2))
        icolor +=1
    plt.legend()
    plt.gca().invert_xaxis()
    #print(f"number of matched sources = {np.sum(zp.matchflag)}")
    
    # add circle
    
    #circle1 = plt.Circle((centerRA, centerDEC), 2*width, color='c',alpha=.2)
    #plt.gca().add_patch(circle1)
    
    plt.savefig(plotprefix+'-se-pan-positions.png')    

def match_coords(secat,pan,verbose=True):
    ###################################
    # match Pan-STARRS1 data to Source Extractor sources
    # remove any objects that are saturated or non-linear in our r-band image
    ###################################


    secoords = SkyCoord(ra=secat['ALPHA_J2000'],dec=secat['DELTA_J2000'],unit=(u.degree,u.degree),frame='icrs')
    #print(pan)
    pancoords = SkyCoord(ra=pan['RAJ2000'],dec=pan['DEJ2000'],unit=(u.degree,u.degree),frame='icrs')
    
    index,dist2d,dist3d = pancoords.match_to_catalog_sky(secoords)
    
    # only keep matches with matched RA and Dec w/in 5 arcsec
    matchflag = dist2d.degree < 5./3600
    
    
    matchedarray1=np.zeros(len(pancoords),dtype=secat.dtype)
    matchedarray1[matchflag] = secat[index[matchflag]]
    
    ###################################################################
    # remove any objects that are saturated, have FLAGS set, galaxies,
    # must have 15 < r < 16.8 for panstarrs mag, based on fitting results
    ###################################################################
    if verbose:
        print(f'\t matched {np.sum(matchflag)} objects')
    fitflag = matchflag  & (pan['rmag'] > 15) & (matchedarray1['FLAGS'] <  1) & (pan['Qual'] < 64)  & (pan['rmag'] < 17.5) #& (matchedarray1['CLASS_STAR'] > 0.95) #& (matchedarray1['MAG_AUTO'] > -11.)
    if verbose:
        print(f'\t number that pass fit {np.sum(fitflag)}')

    return matchedarray1,fitflag
        
def color_correct_panstarrs(pan,filter):
    """
    correcting panstarrs magnitudes into the observed filter systems using conversions from M. Fossati  
    
    Best fit quadratic Intha - PS1_r = 0.0182*(PS1_g-PS1_r)^2 + -0.2662*(PS1_g-PS1_r) + 0.0774
    
    Best fit quadratic Ha4 - PS1_r = 0.0016*(PS1_g-PS1_r)^2 + -0.2134*(PS1_g-PS1_r) + 0.0168

    Best fit quadratic KPSr - PS1_r = 0.0084*(PS1_g-PS1_r)^2 + -0.0420*(PS1_g-PS1_r) + 0.0036

    Best fit quadratic KPHr - PS1_r = 0.0170*(PS1_g-PS1_r)^2 + -0.1864*(PS1_g-PS1_r) + 0.0213

    Best fit quadratic INTSr - PS1_r = 0.0023*(PS1_g-PS1_r)^2 + -0.0122*(PS1_g-PS1_r) + 0.0003


    """
    PS1_r = pan['rmag']
    PS1_g = pan['gmag']
    pan_gr_color = pan['gmag'] - pan['rmag']        
    if filter == 'r':
        print("correcting color for r filter at BOK")                        
        #R = pan['rmag']
        #Best fit quadratic KPSr - PS1_r = 0.0084*(PS1_g-PS1_r)^2 + -0.0420*(PS1_g-PS1_r) + 0.0036
        R = PS1_r + 0.0084*(PS1_g-PS1_r)**2 + -0.0420*(PS1_g-PS1_r) + 0.0036            
        # bok is using the kpno halpha+4nm filter, so use the same correction for these
    elif filter == 'ha':
        print("correcting color for ha filter at KPNO")                        
        #Best fit quadratic Ha4 - PS1_r = 0.0016*(PS1_g-PS1_r)^2 + -0.2134*(PS1_g-PS1_r) + 0.0168
        #R = pan['rmag']
        R = PS1_r + 0.0016*(PS1_g-PS1_r)**2 + -0.2134*(PS1_g-PS1_r) + 0.0168
    return R, pan_gr_color

def get_offsets(secat,pan,filter,verbose=True,naper=5,nsigma=3):
    
    se_matchedarray,fitflag = match_coords(secat,pan,verbose=verbose)
    R,pan_gr_color = color_correct_panstarrs(pan,filter)

    # fit R vs se_matchedarray

    # Start fitting procedure
    flag = fitflag
    bestc = np.array([0,0],'f')
    delta = 100.     
    x = R[flag] # expected mag from panstarrs
    color = pan_gr_color[flag]
    y = se_matchedarray['MAG_APER'][:,naper][flag]
    yerr = se_matchedarray['MAGERR_APER'][:,naper][flag]
    
    #print(f"len(x) = {len(x)}, len(color)= {len(color)}")
    while delta > 1.e-3:
        #c = np.polyfit(x,y,1)
        t = curve_fit(zpfunc,x,y,sigma = yerr)
        # convert to format expected from polyfit
        c = np.array([1.,t[0][0]])
        if verbose:
            print('number of points retained = ',np.sum(flag))
        yfit = np.polyval(c,x)
        residual = (yfit - y)

        delta = abs(bestc[1] - c[1])
        bestc = c
        MAD = mad_std(residual)#1.48*np.median(abs(residual - np.median(residual)))
        #clip_flag = sigma_clip(zim,sigma=3,maxiters=10,masked=True)            
        flag =  (abs(residual - np.median(residual)) < nsigma*MAD)
        if sum(flag) < 2:
            print('WARNING: ONLY ONE DATA POINT LEFT FOR')
            x = x
            y = y
            residual = residual
            sys.exit()
        #flag =  (abs(residual) < nsigma*np.std(residual))
        x = x
        y = y
        residual =residual[flag]
        x = x[flag]
        y = y[flag]
        yerr = yerr[flag]
        color = color[flag]
    zp = c[0]
    yplot = se_matchedarray['MAG_APER'][:,naper][fitflag]
    magfit = np.polyval(bestc,R[fitflag])
    residual_all = 10.**((magfit - yplot)/2.5)        
    residual_all = residual_all
    s = 'residual (mean,std) = %.3f +/- %.3f'%(np.mean(residual_all),np.std(residual_all))
    if verbose:
        print(s)
        
    ###################################
    # Show location of residuals
    ###################################
    plt.figure(figsize=(6,4))
        
    s = ' (mean,std,MAD = {:.2f},{:.2f},{:.2f})'.format(np.mean(residual_all),np.std(residual_all),MAD2(residual_all))
    #s = str(MAD(residual_all))
    plotprefix = 'test'
    plt.title(plotprefix+s)
    residual_allx = se_matchedarray['X_IMAGE'][fitflag]
    residual_ally = se_matchedarray['Y_IMAGE'][fitflag]
    residual_allra = se_matchedarray['ALPHA_J2000'][fitflag]
    residual_alldec = se_matchedarray['DELTA_J2000'][fitflag]
    plt.scatter(se_matchedarray['X_IMAGE'][fitflag],se_matchedarray['Y_IMAGE'][fitflag],c = (residual_all),vmin=.94,vmax=1.04,s=15)
    cb=plt.colorbar()
    cb.set_label('f-meas/f-pan')

    return residual_allx,residual_ally,residual_all,residual_allra,residual_alldec,zp

#################################################################
### MAIN PROGRAM
#################################################################

NAXIS1 = 4032
NAXIS2 = 4096


##
# image name is sent in command line
##
imname = sys.argv[1]


        
ivar_name = imname.replace('ooi','oow').replace('mksb','ksb')

if imname.find('r_v1') > -1:
    image_filter = 'r'
if imname.find('Ha+4nm') > -1:
    image_filter = 'ha'
if imname.find('Ha4nm') > -1:
    image_filter = 'ha'
dq_name = imname.replace('ooi','ood').replace('mksb','ksb')

##
# run SE on MEF image
##
os.system('ln -s ~/github/HalphaImaging/astromatic/default.sex.BOK .')
froot = imname.replace('.fits','')
os.system(f"sex {imname} -c default.sex.BOK  -CATALOG_NAME {froot}.cat")

##
# read in se table
##
secat = froot+'.cat'
shdu = fits.open(secat)

##
# get aprox center RA and DEC
##
centerRA = min(shdu[2].data['ALPHA_J2000'])
centerDEC = min(shdu[2].data['DELTA_J2000'])
width = 2*0.6

##
# get panstarrs catalog
##
pcat = froot+'_pan_tab.csv'
# match each catalog with panstarrs
if os.path.exists(pcat):
    print("found panstarrs cat")
    ptab = Table.read(pcat)
else:
    pan = panstarrs_query(centerRA, centerDEC, width,maxmag=18)
    ptab = Table(pan)
    ptab.write(pcat,format='csv',overwrite=True)


##
# plot positions
# compare SE and panstarrs catalogs
##
if not os.path.exists('plots'):
    os.mkdir('plots')
plot_se_pan_positions(shdu,ptab,plotprefix = 'plots/'+froot)

# get offsets for each amp
se_extensions = [2,4,6,8]
image_extensions = [1,2,3,4]

allx = []
ally = []
allra = []
alldec = []
allresidual = []
allzp = []

for i in se_extensions:
    x,y,residual,ra,dec,zp = get_offsets(shdu[i].data,ptab,image_filter)
    allx.append(x)
    ally.append(y)
    allresidual.append(residual)
    allra.append(ra)
    alldec.append(dec)
    allzp.append(zp)


##
# plot residuals
##

plt.figure(figsize=(10,8))
for i in range(len(allx)):
    plt.scatter(allra[i],alldec[i],c=allresidual[i],vmin=.75,vmax=1.25)
plt.colorbar()
plt.gca().invert_xaxis()
plt.savefig('plots/'+froot+'-residual-xy.png')


##
# Normalize amps
##

##
# get global median for all ccds, so that ccds are normed relative to each other
##

hdu = fits.open(imname) # read in median-subtracted image
ihdu = fits.open(ivar_name) # read in median-subtracted image

allresid1d = [item for sublist in allresidual for item in sublist]
zp_ref = np.mean(np.array(allzp))
global_med = np.median(allresid1d)
print(f"global median for all ccds = {global_med:.3f}")


for h in range(1,len(hdu)):
    ##
    # calculate the offset for each amplifier and scale the data accordingly
    ##
    quad = 0
    print(f"CCD {h}:")
    ccd_med = np.ma.median(allresidual[h-1])
    zp_scale = 10.**((zp_ref - allzp[h-1])/2.5)
    print(f"zp ccd = {allzp[h-1]:.2f}, zpref = {zp_ref:.2f}, zpscale = {zp_scale:.3f}")
    for ix in range(2):
        xmin = 0 + NAXIS1//2*ix
        xmax = NAXIS1//2 * (ix +1)

        for iy in range(2):
            ymin = 0 + NAXIS2//2*iy
            ymax = NAXIS2//2 * (iy+1)
            #print(xmin,xmax,ymin,ymax)
            flag = (allx[h-1] > xmin) & (allx[h-1] < xmax) &\
                (ally[h-1] > ymin) & (ally[h-1] < ymax) 
            amp_med = np.ma.median(allresidual[h-1][flag])

            amp_scale = ccd_med/amp_med * zp_scale

            print(f"\tmedian and scale for quadrant {quad} = {amp_med:.3f} {amp_scale:.3f}")
            ##
            # scale the data so that each ccd/amp has the same ZP
            # scale the image and ivar accordingly
            ##
            hdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*hdu[h].data[ymin:ymax,xmin:xmax]

            ##
            # what is the correct way to scale the weights?  same as image or inverse???
            # in the weight image, high numbers are good
            ##
            ihdu[h].data[ymin:ymax,xmin:xmax] = amp_scale*ihdu[h].data[ymin:ymax,xmin:xmax]

            # add scaling info to the header
            header_key = f"CCD{h}_{quad}"
            header_value = f"{amp_scale:.3f}"
            hdu[0].header.set(header_key,value=float(header_value),comment="PS1 scaling for amp")
            quad += 1

# TODO - should really add scaling factors to the image header so I can compare among images

# save sky-subtracted, scaled image
hdu.writeto('z'+imname,overwrite=True)

ihdu.writeto('zm'+ivar_name,overwrite=True)


shutil.copy(dq_name,'zm'+dq_name)


# check to see if header files exist
# if they do, then copy to have z prefix

scamp_header = 'm'+imname.replace('.fits','.head')
if os.path.exists(scamp_header):
    shutil.copy(scamp_header,'z'+scamp_header)

