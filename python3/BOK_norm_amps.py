#!/usr/bin/env python

import os
import sys

sys.path.append(os.getenv("HOME")+'/github/HalphaImaging/python3/')

import getzp
import ccdproc
import numpy as np

import multiprocessing as mp
from astropy.time import Time
from matplotlib import pyplot as plt
infall_results = []

def process_output(output):
    filename = []
    obstime = []
    amp = []
    zp = []
    zperr = []
    prefix = []
    for t in output:
        # splice image name to get info

        filename.append(t[0])
        if t[0].find('Ha4nm') > -1:
            obstime.append(t[0][14:37])
        else:
            obstime.append(t[0][10:33])
        amp.append(int(t[0].split('_')[1].split('PA')[0]))
        prefix.append(t[0].split('_')[0])
        zp.append(t[1])
        zperr.append(t[2])
    ufiles = set(prefix)
    print(len(ufiles))
    zp = np.array(zp)
    zperr = np.array(zperr)
    amp = np.array(amp)

    scaled_zp = np.ones(len(zp))
    for p in ufiles:
        print(p)
        flag = np.zeros(len(prefix),'bool')
        for i in range(len(prefix)):
            flag[i] = prefix[i] == p

        avezp = np.median(zp[flag])
        scaled_zp[flag] = zp[flag]/avezp
        print(np.sum(flag),avezp)
    return filename,obstime,amp,zp,zperr,prefix,scaled_zp


def make_plots(filename,obstime,amp,zp,zperr,prefix,scaled_zp):
    time = Time(obstime)


    mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    print(len(mycolors))

    delta_time = (time - np.min(time)).value*24

    plt.figure(figsize=(8,6))
    plt.scatter(amp,zp,c=delta_time,alpha=.9,s=20)
    plt.colorbar(label='Time (hr)')
    plt.xlabel('Amplifier',fontsize=16)
    plt.ylabel('ZP (AB)',fontsize=16)
    junk = plt.xticks(np.arange(1,17,2))
    plt.title(f"Images {prefix[0]} -- {prefix[-1]}")
    plt.figure(figsize=(8,6))
    plt.scatter(amp,scaled_zp,c=delta_time,alpha=.9,s=20)
    plt.colorbar(label='Time (hr)')
    plt.xlabel('Amplifier',fontsize=16)
    plt.ylabel('scaled ZP',fontsize=16)
    junk = plt.xticks(np.arange(1,17,2))
    plt.grid()
    plt.ylim(.995,1.005)

    plt.title(f"Images {prefix[0]} -- {prefix[-1]}")    
    chips = set(amp)
    for c in chips:
        flag = amp == c
        print(f"med/ave offset for chip{c:02d} = {np.median(scaled_zp[flag]):.4f},{np.mean(scaled_zp[flag]):.4f}")


def collect_results_zp(result):

    global results
    infall_results.append(result)

def get_images(matchstring="*P.fits"):
    # get images
    ic = ccdproc.ImageFileCollection(os.getcwd(),glob_include=matchstring)
    return ic

def get_central_coord(image):
    # get central coord
    # can get approx from amp 13: add -.15 deg to RA and +.15 to DEC
    h = fits.getheader(image)
    ra = h['CRVAL1']
    dec = h['CRVAL2']    
    return ra-.15,dec+.15

# download panstarrs for full mosaic
# get function from getzp.py - panstarrs_query


def get_zp(image,pfilter,verbose):
    # solve for ZP
    if verbose:
        zpclass,zp,zperr = getzp.main(['--image',image,'--instrument','b','--filter',pfilter,'--verbose'])
    else:
        zpclass,zp,zperr = getzp.main(['--image',image,'--instrument','b','--filter',pfilter])
    return image,zp,zperr


def mp_zp(all_images,verbose=False):
    #######################################################
    # multiprocess this part
    # for each image and each amp, calc zp from panstarrs
    #######################################################
    infall_pool = mp.Pool(mp.cpu_count())
    myresults = [infall_pool.apply_async(get_zp,args=(im,'r',verbose),callback=collect_results_zp) for im in all_images]
    
    infall_pool.close()
    infall_pool.join()
    zp_results = [r.get() for r in myresults]
    return zp_results

    # testing one image at a time to make sure it works ok
    #zp_list=[]
    #zperr_list=[]
    #for im in all_images:
    #    zp,zperr = get_zp(im,'r')
    #    zp_list.append(zp)
    #    zperr_list.append(zperr)
    #return zp_list,zperr_list
def main(matchstring="*P.fits",verbose=False):

    ic = get_images(matchstring=matchstring)
    # print zp for all images, sorted by amplifier

    zpfiles = ic.files_filtered()

    zpresults = mp_zp(zpfiles,verbose=verbose)

    return zpresults
    # normalize all images to the same ZP


    # subtract sky???
