#!/usr/bin/env python
"""

USAGE

just call it in the directory containing the ksb images

"""

import os
import subprocess
import sys
import glob
import time
import multiprocessing as mp
HOME = os.getenv("HOME")

infall_results = []
def collect_results(result):

    global results
    infall_results.append(result)

def processone(image):
    #print(f"processing {i}/{len(filelist)}")
    
    
    #image = filelist[i]
    outimage = 'z'+image
    #print("processing ",image)
    if os.path.exists(outimage):
        print(f'{outimage} found - moving to next object')
        return
    else:
        print(f"processing {image} -> {outimage}")

    # updating to use the new program
    program= f"{HOME}/github/HalphaImaging/python3/BOK_fixampoffsets.py"    
    cmd = f"python {program} {image}"
    os.system(cmd)
    return 1
    #cmds = ['python3', program,image]
    #print(f"Submitting job to process {input_file}")
    #process = subprocess.Popen(cmds)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #stdout, stderr = process.communicate()
    #print(stdout.decode())
    #print(stderr.decode())
def processall():

    # changing prefix from ksb to mksb b/c this is run on median-subtracted images
    filelist = glob.glob('mksb*ooi*_v1.fits')
    filelist.sort()
    print(f'found {len(filelist)} files to process...')
    #t = [image1 for image1 in filelist]
    #print(t)
    # set up multiprocessing pool
    image_pool = mp.Pool(mp.cpu_count())
    myresults = [image_pool.apply_async(processone,args=(im,),callback=collect_results) for im in filelist]
    
    image_pool.close()
    image_pool.join()
    infall_results = [r.get() for r in myresults]
    return infall_results

if __name__ == '__main__':
    infall_results = processall()
