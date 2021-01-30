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

def subtract_median_sky(data):
    ''' subtract median sky from image data '''
    mask = make_source_mask(data,nsigma=3,npixels=5,dilate_size=10)
    mean,median,std = sigma_clipped_stats(data,sigma=3.0,mask=mask)
    data -= median
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
        newheader.set(keys[i],value=values[i],comment=comments[i])

    return newheader
                            
            
        
        

