#!/usr/bin/env python

"""
GOAL:
- generate cutouts for images with halpha cutouts

NOTES:
- run in the cutouts directory

- it will go into each subdir and run the script to generate cutouts!

TODO : add a random wait time before starting program so that image servers are not overwhelmed

"""
#import glob
import os
import sys
import time
import np
homedir = os.getenv('HOME')
sys.path.append(homedir+'/github/HalphaImaging/python3/')

randomtimes = 
def getone(d):
    import glob    
    workingdir = os.getcwd()
    os.chdir(d)

    imlist = glob.glob('VF*-R.fits')         
    #imlist = a + b + c + d + e
    #imlist = glob.glob(d+'*-R.fits')
    if len(imlist) == 0:
        print("did not find any images!")
        return
    else:
        print(f"found {len(imlist)} r-band images")
        
    for im in imlist:
        #####################################################################################
        # add a few seconds sleep time to avoid overwhelming servers when running in parallel
        #####################################################################################        
        sleeptime = 5.*np.random.random()
        time.sleep(sleeptime)

        #####################################################################################
        # call program to download images
        #####################################################################################        
        print('python '+homedir+'/github/HalphaImaging/python3/plot_cutouts_ha.py --r '+im+' --plotall')
        os.system('python '+homedir+'/github/HalphaImaging/python3/plot_cutouts_ha.py --r '+im+' --plotall')
        #####################################################################################        
    os.chdir(workingdir)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description ='get cutout images for all cutouts.  Run this from the directory that contains subdirectories for each galaxy.  For example: /data-pool/Halpha/halphagui-output-20230711/')
    parser.add_argument('--onegal',dest = 'onegal', default=None, help = "use this to run on one image only.  Specify the full path to the r-band image. ")
    args = parser.parse_args()



    if args.onegal is not None:
        # run on one image only
        # this is used for running in parallel
        getone(args.onegal)
    else:
        import multiprocessing as mp
        dirlist = glob.glob('VFID*')
        dirlist.sort()
        workingdir = os.getcwd()
        dirlist = [f.path for f in os.scandir(workingdir) if f.is_dir() ]
        

    
        for d in dirlist:
            image_pool = mp.Pool(mp.cpu_count())
            #image_pool = mp.pool.ThreadPool(mp.cpu_count())    
            myresults = [image_pool.apply_async(getone,args=(d,)) for d in dirlist]
            
            image_pool.close()
            image_pool.join()
            image_results = [r.get() for r in myresults]
    
        

