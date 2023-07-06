#!/usr/bin/env python
"""
GOAL: run getzp, filter ratio, build psf, and build html

USAGE:

python process_one_coadd.py rimage haimage

or 

python ~/github/HalphaImaging/python3/process_one_coadd.py VF-177.200+56.055-INT-20220502-VFID0957-r-shifted.fits  VF-177.200+56.055-INT-20220502-VFID0957-Halpha.fits


start in the directory containing the images

"""

import sys
import os

rimage = sys.argv[1]
haimage = sys.argv[2]

startdir = os.getcwd()

if 'BOK' in rimage:
    instrument='b'
elif 'INT' in rimage:
    instrument='i'
elif 'HDI' in rimage:
    instrument='h'
elif 'MOS' in rimage:
    instrument='m'

# Solve for zp
#s = f"python ~/github/HalphaImaging/python3/getzp.py --image {rimage} --instrument {instrument} --filter r --flatten 1"
# running for the rogue INT images in the all-virgo-cutouts, so don't need to flatten again
if instrument == 'm':
    s = f"python ~/github/HalphaImaging/python3/getzp.py --image {rimage} --instrument {instrument} --filter R"
else:
    s = f"python ~/github/HalphaImaging/python3/getzp.py --image {rimage} --instrument {instrument} --filter r"
os.system(s)

#s = f"python ~/github/HalphaImaging/python3/getzp.py --image {haimage} --instrument {instrument} --filter ha --flatten 1"
s = f"python ~/github/HalphaImaging/python3/getzp.py --image {haimage} --instrument {instrument} --filter ha "
os.system(s)


if not os.path.exists('SEcats_getzp'):
    os.mkdir('SEcats_getzp')

images = [rimage,haimage]
for im in images:
    se_cat = im.replace('.fits','.cat')

    s = f"mv {se_cat} SEcats_getzp/."
    os.system(s)

s = f"python ~/github/halphagui/batch_filterratio.py --oneimage {rimage}"
os.system(s)


# build psf images

os.chdir('/data-pool/Halpha/psf-images/')

for im in images:
    full_path_rimage = os.path.join(startdir,im)
    s = f"python ~/github/halphagui/buildpsf.py --image {full_path_rimage} "
    if instrument == 'i':
        s += " --int"
    elif instrument == 'b':
        s += " --bok"
    os.system(s)


# build the webpage
os.chdir('/data-pool/Halpha/html_dev/')
full_path_rimage = os.path.join(startdir,rimage)
s = f"python ~/github/Virgo/programs/build_web_coadds2.py --oneimage {full_path_rimage} "
os.system(s)
