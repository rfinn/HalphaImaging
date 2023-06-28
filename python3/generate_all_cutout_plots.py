#!/usr/bin/env python

"""
GOAL:
- generate cutouts for images with halpha cutouts

NOTES:
- run in the cutouts directory

- it will go into each subdir and run the script to generate cutouts!

"""
#import glob
import os
import multiprocessing as mp
import sys
homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/HalphaImaging/python3/')

dirlist = glob.glob('VFID*')
dirlist.sort()
workingdir = os.getcwd()
dirlist = [f.path for f in os.scandir(workingdir) if f.is_dir() ]
def getone(d):
    imlist = glob.glob(d+'*-R.fits')
    if len(imlist) == 0:
        imlist = glob.glob(d+'*-r.fits')
    if len(imlist) == 0:
        continue
        
    for im in imlist:
        os.system('python '+homedir+'/github/HalphaImaging/python3/plot_cutouts_ha.py --r '+im+' --plotall')
    os.chdir(workingdir)    
for d in dirlist:
    image_pool = mp.Pool(mp.cpu_count())
    #image_pool = mp.pool.ThreadPool(mp.cpu_count())    
    myresults = [image_pool.apply_async(buildone,args=(d,)) for d in dirlist]
    
    image_pool.close()
    image_pool.join()
    image_results = [r.get() for r in myresults]
    
        

