#!/usr/bin/python
# coding: utf-8

import os
import sys
import glob
import argparse
import numpy as np
from astropy.io import fits
from scipy.stats import mode
import matplotlib.pyplot as plt


parser = argparse.ArgumentParser(description ='Run sextractor on mosaics')
parser.add_argument('--d',dest = 'd', default =' ~/github/HalphaImaging/astromatic', help = 'Locates path of config files')
parser.add_argument('--c',dest = 'c', default ='all', help = 'Cluster name to find scale of. Defaults to all in directory')
parser.add_argument('--ns',dest = 'ns', default =True, action = 'store_false', help = 'If included will not run sextractor')

args = parser.parse_args()
os.system('cp ' +args.d + '/default.* .')

def FindScale(cluster):
    print "Cluster Name:", cluster
    hacat = glob.glob(cluster + "_ha*.cat")
    rcat = glob.glob(cluster + "_R.cat")
    try:
        hacat = hacat[0]
        rcat = rcat[0]
    except IndexError:
        print "No catalogs found"
        print
        return
    hadat = fits.getdata(hacat,2)
    rdat = fits.getdata(rcat,2)
    haflux = hadat["FLUX_AUTO"]
    rflux = rdat["FLUX_AUTO"]
    hflags = hadat["FLAGS"]
    rflags = rdat["FLAGS"]
    keepflag = np.ones(len(rflags),dtype = bool)
    for i in range(len(hflags)):
        if (rflags[i]+hflags[i])>0: # Doesn't accept any flags
            keepflag[i] = False
    rgflx = rflux[keepflag]
    hagflx = haflux[keepflag]
    qflux = hagflx/rgflx
    #nflag = np.ones(len(qflux),dtype = bool)
    nflag = (qflux < .15)
    #for i in range(len(qflux)):
    #    if abs(qflux[i]) >.15: # Defines what is an outlier
    #        nflag[i] = False
    qflux = qflux[nflag]
    plt.figure()
    plt.hist(qflux,bins = np.arange(.03,.055,.002))
    plt.title("Halpha Flux / R Flux")
    plt.xlabel("Scale")
    plt.ylabel("Amount")
    plt.savefig(cluster+'_scale.png')
    print "Number of Points Found:", len(rflags)
    print "Number of Uncompromised Points:", len(rgflx)
    print "Number of Points w/o Outliers:", len(qflux)
    print "Average Scale:", np.mean(qflux)
    print "Median Scale:", np.median(qflux)
    print 
    return np.mean(qflux)



if args.c == 'all':
    if args.ns:
        hafiles = glob.glob('*_ha*')
        hafiles = set(hafiles) - set(glob.glob("*.*"))
        for it in hafiles:
            t = it.split('_')
            print "Running SExtractor on",t[0]
            ir = t[0]+'_R'
            os.system('/usr/bin/sextractor ' + it + '.coadd.fits -c default.sex.hdi -CATALOG_NAME ' + it + '.cat')
            os.system('/usr/bin/sextractor ' + it + '.coadd.fits,' + ir + '.coadd.fits -c default.sex.hdi -CATALOG_NAME ' + ir + '.cat')
            print
            print
    
    clusters = glob.glob('*_ha*')
    clusters = set(clusters) - set(glob.glob("*.*"))
    for i in clusters:
        t = i.split('_')
        FindScale(t[0])
else:
    if args.ns:
        cluster = glob.glob(args.c+'_ha*')
        try:
            cl = cluster[0]
        except IndexError:
            print "No catalogs found for cluster", args.c
            sys.exit(0)
        
        cn = cl.split('_')
        print 'Running SExtractor on',cn[0]
        cr = cn[0]+'_R'
        os.system('/usr/bin/sextractor ' + cl + '.coadd.fits -c default.sex.hdi -CATALOG_NAME ' + cl + '.cat')
        os.system('/usr/bin/sextractor ' + cl + '.coadd.fits,' + cr + '.coadd.fits -c default.sex.hdi -CATALOG_NAME ' + cr + '.cat')
        print
        print
    
    FindScale(args.c)
        



