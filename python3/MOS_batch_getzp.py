#!/usr/bin/env python
"""
Run from the directory containing the VF*R.fits images
"""
import glob
import os
from matploblib import pyplot as plt

Rfiles = glob.glob('VF*R.fits')
Rfiles.sort()
for rf in Rfiles:
    hf = rf.replace('R.fits','Ha4.fits')
    str1 = f"python ~/github/HalphaImaging/python3/getzp.py --image {rf} --instrument m --filter R --flatten 1 --spline"
    print(str1)
    os.system(str1)
    plt.close("all")
    str2 = f"python ~/github/HalphaImaging/python3/getzp.py --image {hf} --instrument m --filter ha --flatten 1 --spline"
    print(str2)
    os.system(str2)
    plt.close("all")




    
    
