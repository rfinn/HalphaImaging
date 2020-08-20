#!/usr/bin/env python

'''
Read in sextractor catalog and plot image with detected objects

'''

import numpy as np
from matplotlib import pyplot as plt
import argparse

parser = argparse.ArgumentParser(description ='Get cutouts for NSA galaxies within field of view of mosaic and redshift range of designated H-alpha filter')
parser.add_argument('--image', dest = 'image', default = 'halpha.fits', help = 'mosaic/HDI image to make cutouts from')
