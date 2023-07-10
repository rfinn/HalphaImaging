#!/usr/bin/env python

'''
https://photutils.readthedocs.io/en/stable/background.html


REFERENCES:
masking
https://mwcraig.github.io/ccd-as-book/08-02-Creating-a-mask.html

'''

from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import ccdproc
from photutils import make_source_mask
from astropy.io.fits import Header
import numpy as np

def subtract_median_sky(data,getstd=False,getmedian=True,subtract=True):
    ''' subtract median sky from image data '''
    mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=5)
    masked_data = np.ma.array(data,mask=mask)
    #clipped_array = sigma_clip(masked_data,cenfunc=np.ma.mean)

    mean,median,std = sigma_clipped_stats(masked_data,sigma=3.0,cenfunc=np.ma.mean)
    if subtract:
        data -= median
    if getstd:
        return data,median,std
    
    elif getmedian:
        return data,median
    
    else:
        return data
    
def make_ccdmask(flat1,flat2=None):
    ''' make bad pixel mask from flat image.  use ratio of flats if flat2 is given '''
    if flat2 is None:
        flatimage = fits.getdata(flat1)
    else:
        flat1 = fits.getdata(flat1)
        flat2 = fits.getdata(flat2)
        flatimage = flat1/flat2

    # use ccdmask to identify bad pixels

    t = ccdproc.ccdmask(flatimage)
    return t

def convert_headfile_header(filename):
    f = open(filename,'r')
    newheader = Header()
    keys=[]
    values=[]
    comments=[]
    for line in f:
        if line.startswith('COMMENT'):
            continue
        if line.find('=') > -1:
            t = line.split('=')
            keys.append(t[0].strip())
            m = t[1].split('/')
            try:
                values.append(float(m[0]))
            except ValueError:
                values.append(m[0].strip().replace("'",""))
    
            comments.append(m[1].strip())
    f.close()
    # construct new header
    for i in range(len(keys)):
        try:
            newheader.set(keys[i],value=values[i],comment=comments[i])
        except ValueError:
            # getting error ValueError: Floating point inf values are not allowed in FITS headers.
            print('WARNING: did not like header entry: ',keys[i],values[i])
            continue
    return newheader
                            
def fit_surface2points(x,y,z):
    """would like another way to fit residuals, maybe boxcar smooth?
    polynomial sometimes introduces too much structure

    https://stackoverflow.com/questions/39727040/matplotlib-2d-plot-from-x-y-z-values
    """
    from scipy.interpolate import interp2d

    # f will be a function with two arguments (x and y coordinates),
    # but those can be array_like structures too, in which case the
    # result will be a matrix representing the values in the grid 
    # specified by those arguments
    f = interp2d(x,y,z,kind="linear")

    x_coords = np.arange(min(x),max(x)+1)
    y_coords = np.arange(min(y),max(y)+1)
    Z = f(x_coords,y_coords)

    fig = plt.imshow(Z,
           extent=[min(x_list),max(x_list),min(y_list),max(y_list)],
           origin="lower")
    # Show the positions of the sample points, just to have some reference
    fig.axes.set_autoscale_on(False)
    plt.scatter(x_list,y_list,400,facecolors='none')
    pass    
        
        

