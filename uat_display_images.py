#!/usr/bin/env python




import glob
import argparse

from matplotlib import pyplot as plt
import numpy as np
from scipy.stats import scoreatpercentile
from astropy.io import fits

parser = argparse.ArgumentParser(description ='Groups images by filter and creates flatfield images')
parser.add_argument('--filestring', dest='filestring', default='ztr', help='match string for input files (default =  ztr, which gets ztr*.fits)')
#parser.add_argument('--', dest='pixelscalex', default='0.00011808', help='pixel scale in x (default = 0.00011808)')
args = parser.parse_args()
files = sorted(glob.glob(args.filestring+'*.fits'))
nfiles=len(files)

for f in files:
    plt.close('all')
    plt.figure()
    imdat = fits.getdata(f)
    print('Displaying image '+f)
    vmin,vmax=scoreatpercentile(imdat,[.5,99.9])
    plt.imshow(imdat,origin='upper',vmin=vmin, vmax=vmax, cmap = 'Greys')
    plt.show()
    t = raw_input('press any key for next image; press q to quit \n')
    if t.find('q') > -1:
        break
