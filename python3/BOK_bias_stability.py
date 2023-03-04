#!/usr/bin/env python

"""
GOAL: test bias levels to see if they are stable with time

PROCEDURE:
* get list of bias images
* for each image, get time of exposure
* for each hdu, get median, mean, std in counts
* save counts in an Nimage x 16 array
* plots counts in each amplifier vs time

* Run in directory containing raw bias images

"""
import ccdproc
import os
from astropy.io import fits
from astropy.time import Time
import numpy as np

from matplotlib import pyplot as plt



ic = ccdproc.ImageFileCollection(os.getcwd(), glob_include='*.fit*')



bias_files = ic.files_filtered(imagetyp="zero")
print('Parsing {} files'.format(len(bias_files)))
# number of bias frames
nbias = len(bias_files)

# get number of extensions
t = fits.open(bias_files[0])
nextension = len(t) -1
t.close()

# make data frame to hold bias of each amp in each image
bias_mean = np.zeros((nbias,nextension))
bias_median = np.zeros((nbias,nextension))
bias_std = np.zeros((nbias,nextension))
time = [] # list of times

for i,bias in enumerate(bias_files):
    #print(i,bias)
    hdu = fits.open(bias)
    # record the time of the exposure
    test_time = hdu[0].header['TIME-OBS']
    test_date = hdu[0].header['DATE-OBS']
    newtime = test_date+" "+test_time    
    time.append(newtime)
    # skip the primary hdu, which doesn't have an associated image
    for j in range(1,len(hdu)):
        #print(i,j)
        bias_mean[i,j-1] = np.mean(hdu[j].data)
        bias_median[i,j-1] = np.median(hdu[j].data)
        bias_std[i,j-1] = np.std(hdu[j].data)
    hdu.close()

time = Time(time)


mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
print(len(mycolors))
delta_time = (time - time[0]).value*24

# sort data points by time
isorted = np.argsort(delta_time)

plt.figure(figsize=(10,6))


for i in range(nextension):
    if i >= len(mycolors):
        myls='--'
    else:
        myls='-'
    plt.subplot(1,2,1)
    plt.plot(delta_time[isorted],bias_mean[:,i][isorted],label=f"ext={i}",marker='*',ls=myls)
    plt.subplot(1,2,2)
    plt.plot(delta_time[isorted],bias_median[:,i][isorted],label=f"ext={i}",marker='^',ls=myls)

plt.subplot(1,2,1)
plt.title("Mean Counts",fontsize=18)
#plt.legend()
plt.grid()
plt.ylabel("Counts (ADU)",fontsize=16)
plt.ylim(850,1175)
plt.xlabel("Time (hrs)",fontsize=16)

plt.subplot(1,2,2)
plt.legend(bbox_to_anchor=(1.05, 1),ncol=1)
plt.grid()
plt.xlabel("Time (hrs)",fontsize=16)
plt.title("Median Counts",fontsize=18)
plt.ylim(850,1175)
plt.savefig('bias_vs_time.png')


mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
print(len(mycolors))
delta_time = (time - np.min(time)).value*24

# sort data points by time
isorted = np.argsort(delta_time)

plt.figure(figsize=(10,6))


for i in range(nextension):
    if i >= len(mycolors):
        myls='--'
    else:
        myls='-'
    plt.subplot(1,2,1)
    plt.plot(delta_time[isorted],bias_mean[:,i][isorted]/bias_mean[:,i][isorted][0],label=f"ext={i}",marker='*',ls=myls)
    plt.subplot(1,2,2)
    plt.plot(delta_time[isorted],bias_median[:,i][isorted]/bias_median[:,i][isorted][0],label=f"ext={i}",marker='^',ls=myls)

plt.subplot(1,2,1)
plt.title("Mean Counts",fontsize=18)
#plt.legend()
plt.grid()
plt.ylabel("Counts/Counts(t0)",fontsize=16)
#plt.ylim(850,1175)
plt.xlabel("Time (hrs)",fontsize=16)

plt.subplot(1,2,2)
plt.legend(bbox_to_anchor=(1.05, 1),ncol=1)
plt.grid()
plt.xlabel("Time (hrs)",fontsize=16)
plt.title("Median Counts",fontsize=18)
#plt.ylim(850,1175)
plt.savefig('bias_vs_time_normed.png')
