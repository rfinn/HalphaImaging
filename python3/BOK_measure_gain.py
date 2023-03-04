#!/usr/bin/env python

"""
GOAL: measure the gain in each amplifier

PROCEDURE:
* get list of flat images
* for each image, get time of exposure
* for each hdu, get median, mean, std in counts in 100 random locations
* plots noise vs counts in each amplifier
* plot noise/counts vs time in each amplifier to see if there is drift

* Run in directory containing bias-subtracted flats

"""
import ccdproc
import os
from astropy.io import fits
from astropy.time import Time
import numpy as np

from matplotlib import pyplot as plt

import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--comb', action='store_true', help='Run gain on combined images.  This assumes uncombined images are *P.fits, and combined images do not have P (theli naming convention).')
parser.add_argument('--ncomb', dest='ncomb',default=1, help='number of images that were combined to make the flat image.  Default is one.')
args = parser.parse_args()

if args.comb:

    ic = ccdproc.ImageFileCollection(os.getcwd(), glob_include="*.fits",glob_exclude='*P.fits')
else:
    ic = ccdproc.ImageFileCollection(os.getcwd(),glob_include='*P.fits')

# once theli processes the flats, they are split into individual amplifiers
flat_files = ic.files_filtered()# imagetyp is not in processed theli images imagetyp="flat")

# cut for testing
#flat_files = flat_files[0:64]

print('Parsing {} files'.format(len(flat_files)))
# number of bias frames
nflat = len(flat_files)

# get number of extensions
nextension = 16
nsamples = 100
chipbuffer = 100
nrootimages = int(nflat/nextension)
# square for measuring stats
regsize = 25
# make data frame to hold bias of each amp in each image
fmean = np.zeros((nrootimages,nextension,nsamples))
fmedian = np.zeros((nrootimages,nextension,nsamples))
fstd = np.zeros((nrootimages,nextension,nsamples))
ftime = [] # list of times

# need to keep track of root image name, each having 16 associate files
imroot = []

for i,mfile in enumerate(flat_files):
    hdu = fits.open(mfile)
    ext = mfile.split('_')
    if ext[0] not in imroot:
        imroot.append(ext[0])
        print(f"processing extensions from {ext[0]}")
        xmax,ymax = hdu[0].data.shape
        
        #print(i,mfile,ext)

        # record the time of the exposure
        if args.comb:
            newtime = '2022-04-23 12:00:00.0'
        else:
            newtime = hdu[0].header['DATE-OBS']

        ftime.append(newtime)
        
        if len(imroot) == 1:
            nimage = 0
        else:
            nimage += 1
    if args.comb:
        ext = int(ext[1].split('.fits')[0])
    else:
        ext = int(ext[1].split('P')[0])


    xmax,ymax = hdu[0].data.shape
    if nimage == 0:
        print(f"ext{ext}: size of ccd is {xmax}x{ymax}")
    xpos = np.random.randint(low=chipbuffer,high=xmax-chipbuffer,size=nsamples)
    ypos = np.random.randint(low=chipbuffer,high=ymax-chipbuffer,size=nsamples)        
    for k in range(len(xpos)):
        #print(f"{xpos-regsize}:{xpos+regsize},{ypos-regsize}:{ypos+regsize}")
        x1 = xpos[k]-regsize
        x2 = xpos[k]+regsize
        y1 = ypos[k]-regsize
        y2 = ypos[k]+regsize
        fmean[nimage,ext-1,k] = np.mean(hdu[0].data[x1:x2,y1:y2])
        fmedian[nimage,ext-1,k] = np.median(hdu[0].data[x1:x2,y1:y2])
        fstd[nimage,ext-1,k] = np.std(hdu[0].data[x1:x2,y1:y2])
        #print(i,ext,k,fmean[i,ext-1,k])
    hdu.close()

time = Time(ftime)


mycolors = plt.rcParams['axes.prop_cycle'].by_key()['color']
print(len(mycolors))

delta_time = (time - np.min(time)).value*24


# sort data points by time
#isorted = np.argsort(delta_time)

#delta_time = delta_time[isorted]
#fmean = fmean[isorted]
#fmedian = fmedian[isorted]
#fstd = fstd[isorted]

def plot_gain(counts,std,time,title=None,ymin=None,ymax=None):
    ax1=plt.figure(figsize=(10,10))
    nflat, nextension, nsamples = counts.shape

    allgains = []
    allax = []    
    for i in range(nextension):
        if i >= len(mycolors):
            mymarker='^'
        else:
            mymarker='^'
        plt.subplot(4,4,i+1)
        
        # plot noise^2 vs counts
        for j in range(nflat):
            plt.scatter(counts[j,i,:],std[j,i,:]**2,c=time[j]*np.ones(nsamples),label=f"ext={i}",alpha=.2,marker=mymarker,vmin=np.min(time),vmax=np.max(time))
        allax.append(plt.gca())
        # calculate the gain
        allgains.append([std[:,i,:].flatten()**2/counts[:,i,:].flatten()])
        
        if i == 0 and title is not None:
            plt.title(f"{title}",fontsize=18)
    
        if i % 4 == 0:
            plt.ylabel("$\sigma^2 \ (ADU)$")
        if i > 11:
            plt.xlabel("Counts (ADU)",fontsize=16)
        plt.gca().set_yscale('log')
        if ymin is not None:
            plt.ylim(ymin,ymax)

        plt.text(.2,.9,f"ext={i+1}",transform=plt.gca().transAxes,horizontalalignment='left')
    plt.colorbar(ax=allax,fraction=.08,label='Time')
    # plot histograms
def plot_gain_hists(counts,std,time,title=None,ymin=None,ymax=None):
    nflat, nextension, nsamples = counts.shape

    allgains = []
    allax = []
    gain = std**2/counts
    for i in range(nextension):
        allgains.append([gain[:,i,:].flatten()])    
    plt.figure(figsize=(10,10))    
    mybins = np.linspace(.5,2,100)
    print('got here. now making hists...')
    for i in range(len(allgains)):
        #print(f'working on panel {i}')
        plt.subplot(4,4,i+1)
        if i%4 != 0:
            plt.yticks([],[])
        if i < 11:
            plt.xticks([],[])

        plt.hist(allgains[i],bins=mybins)
        plt.axvline(x=np.median(allgains[i]),ls='--',color='k',label='med')
        plt.axvline(x=np.mean(allgains[i]),ls='-',color='0.5',label='mean')        
        plt.text(.05,.9,f"med={np.median(allgains[i]):.2f}",transform=plt.gca().transAxes,horizontalalignment='left')
        plt.legend(loc='lower left')
        plt.title(f"amp-{i+1}")
        plt.xlim(min(mybins),max(mybins))
def plot_time_figs():
    plt.figure(figsize=(10,6))


    for i in range(nextension):
        if i >= len(mycolors):
            myls='--'
        else:
            myls='-'
        plt.subplot(1,2,1)
        plt.plot(delta_time[isorted],fmean[:,i][isorted],label=f"ext={i}",marker='*',ls=myls)
        plt.subplot(1,2,2)
        plt.plot(delta_time[isorted],fmedian[:,i][isorted],label=f"ext={i}",marker='^',ls=myls)

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
    plt.savefig('flat_vs_time.png')



    plt.figure(figsize=(10,6))


    for i in range(nextension):
        if i >= len(mycolors):
            myls='--'
        else:
            myls='-'
        plt.subplot(1,2,1)
        plt.plot(delta_time[isorted],fmean[:,i][isorted]/fmean[:,i][isorted][0],label=f"ext={i}",marker='*',ls=myls)
        plt.subplot(1,2,2)
        plt.plot(delta_time[isorted],fmedian[:,i][isorted]/fmedian[:,i][isorted][0],label=f"ext={i}",marker='^',ls=myls)

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
    plt.savefig('flat_vs_time_normed.png')
