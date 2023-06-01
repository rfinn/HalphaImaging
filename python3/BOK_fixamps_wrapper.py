#!/usr/bin/env python

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
    outimage = 'zm'+image
    print("processing ",image)
    #if os.path.exists(outimage):
    #    print(f'{outimage} found - moving to next object')
    #else:
    #    print(f"processing {image} -> {outimage}")
    
    program= f"{HOME}/github/HalphaImaging/python3/BOK_pipeline_fixampoffsets.py"    
    cmd = f"python {program} {image} r"
    os.system(cmd)
    return 1
    #cmds = ['python3', program,image]
    #print(f"Submitting job to process {input_file}")
    #process = subprocess.Popen(cmds)#, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #stdout, stderr = process.communicate()
    #print(stdout.decode())
    #print(stderr.decode())
def processall():

    filelist = glob.glob('ksb*ooi*.fits')
    filelist.sort()
    print(f'found {len(filelist)} files to process...')

    # set up multiprocessing pool
    image_pool = mp.Pool(mp.cpu_count())
    myresults = [image_pool.apply_async(processone,args=(image1),callback=collect_results) for image1 in filelist]
    
    image_pool.close()
    image_pool.join()
    infall_results = [r.get() for r in myresults]

if __name__ == '__main__':
    processall()
