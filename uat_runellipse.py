#!/usr/bin/env python
'''
DESCRIPTION:

Run iraf.stsdas.analysis.isophote.ellipse on R-band image.  Take the parameters
from the R-band fit and use this as input for running ellipse on H-alpha image.



PROCEDURE:

- pass in r_image and Ha_image
- assumes mask is named e.g. CLUSTER-NSAID-mask.fits, e.g. A1367-140231-mask.fits 


USEAGE:

run in ipython -pylab

   %run ~/github/HalphaImaging/uat_runellipse.py --r r_image.fits --ha ha_image.fits

'''

import sys
import os
from matplotlib import pyplot as plt 
from pyraf import iraf
from astropy.io import fits
import ds9
import argparse


parser = argparse.ArgumentParser(description ='Run ellipse on R and Halpha images')
parser.add_argument('--cluster',dest = 'cluster', default = None, help = 'Cluster prefix of image to run ellipse on.')
parser.add_argument('--id',dest = 'id', default = None, help = 'NSAID of image to run ellipse on.')

#parser.add_argument('--r', dest = 'r', default = None, help = 'R-band image')
#parser.add_argument('--ha', dest = 'ha', default = None, help = 'Halpha image')
#parser.add_argument('--scale', dest = 'scale', default = 0.06, help = 'R-band image')
args = parser.parse_args()


#if args.threshold:
#    print 'got a threshold value',args.threshold
#    deblend = float(args.threshold)
#if args.catalog:
#    print 'got an input default.sex file = ',args.catalog
#    defaultcat=args.catalog



#if args.display:
#    import ds9

def call_measure_disk(input_image,mask_image=None,ipa=0,xcenter=None,ycenter=None,minr=2,initialr=15,iellip = .05,wave_band=0,keepfixed=0):
    print 'mask image = ',mask_image
    fdulist = fits.open(input_image)
    t=fdulist[0].data
    n2,n1=t.shape
    if xcenter == None:
        xc=n1/2.
        yc=n2/2.
    else:
        xc=xcenter
        yc=ycenter
    fdulist.close()

    if ipa < -90:
        ipa=ipa+180
    elif ipa > 90:
        ipa = ipa-180
        
    if initialr < 5:
        initialr = 5
    if initialr > 100:
        initialr = 50

    maxr=min(n1,n2)/2.
    if initialr > maxr:
        initialr = 0.25*maxr
    if wave_band == 1:
        initialr=20
    #if maxr > maxrad:
        #maxr=maxrad-5
    # make up a value for now
    mag_zp=25.
    try:
        d.set('frame delete all')
    except:
        d=ds9.ds9()

    d.set('frame 1')
    d.set('single')
    # display r-band image if running Halpha image
    if wave_band == 1:
        d.set('frame 2')
        s='file '+input_image
        try:
            d.set(s)
        except:
            print 'trouble loading ',input_image

    print xc,yc,iellip,ipa,initialr,minr,maxr,mag_zp
    if keepfixed:
        newxcenter,newycenter,newellip,newPA=measure_disk(input_image,xc,yc,ipa,iellip,initialr,minr,maxr,mag_zp,band=wave_band,keepfixed=1,mask_image=mask_image)
        t=raw_input('enter x to quit, i to rerun in interactive mode, any other key to continue\n')
        if t.find('x') > -1:
            return 1
        if t.find('i') > -1:
            newxcenter,newycenter,newellip,newPA=measure_disk(i,disk_image,newxcenter,newycenter,newPA,newellip,initialr,minr,maxr,mag_zp,band=wave_band,keepfixed=1,mask_image=mask_image,einteractive=1)

    else:

        newxcenter,newycenter,newellip,newPA=measure_disk(input_image,xc,yc,ipa,iellip,initialr,minr,maxr,mag_zp,band=wave_band,keepfixed=0,mask_image=mask_image)
        t=raw_input('enter x to quit, i to rerun in interactive mode, any other key to continue\n')
        if t.find('x') > -1:
            return 1
        if t.find('i') > -1:
            newxcenter,newycenter,newellip,newPA=measure_disk(input_image,newxcenter,newycenter,newPA,newellip,initialr,minr,maxr,mag_zp,band=wave_band,keepfixed=0,einteractive=1,mask_image=mask_image)
        newxcenter,newycenter,newellip,newPA=measure_disk(input_image,newxcenter,newycenter,newPA,newellip,initialr,minr,maxr,mag_zp,band=wave_band,keepfixed=1,einteractive=0,mask_image=mask_image)

    return newxcenter,newycenter,newellip,newPA

def measure_disk(input_image,xcenter,ycenter,ipa,iellip,initialr,minr,maxr,zp,band=0,nframe=1,myradius=20,keepfixed=0,einteractive=0,mask_image=None):
    print 'inside measure disk'
    print 'mask image = ',mask_image
    #myradius=initialr
    recentervalue='yes'
    interactivevalue='yes'

    print '\n Running ellipse on ',input_image,'\n'
    print ' Enter ? for a list of options \n'
    print ' Enter c to continue fitting ellipses automatically \n'
    mfile=input_image
    
    xcenter_cutout=xcenter
    ycenter_cutout=ycenter
    
    t=mfile.split('.')
    efile=t[0]+'.tab'
    efile_ascii=t[0]+'.dat'
    imprefix=t[0]
    print 'Running ellipse to fit isophotes to galaxy:',mfile
    if os.path.exists(efile):
        os.remove(efile)
    if keepfixed:
        print input_image,efile,xcenter_cutout,ycenter_cutout,initialr,minr,maxr,ipa
        iraf.ellipse(input=input_image,output=efile,dqf=mask_image,x0=xcenter_cutout,y0=ycenter_cutout,hcenter='yes',recenter=recentervalue,sma0=initialr,minsma=minr,maxsma=maxr,pa=ipa,hpa='yes',ellip=iellip,hellip='yes',interactive='yes',step=0.1,linear='no')
                
    else:
        print "First pass, letting PA and e vary"
        print input_image,efile,xcenter_cutout,ycenter_cutout,initialr,minr,maxr,ipa,
        if mask_image == None:
            iraf.ellipse(input=input_image,output=efile,x0=xcenter_cutout,y0=ycenter_cutout,hcenter='no',recenter='yes',sma0=initialr,minsma=minr,maxsma=maxr,pa=ipa,hpa='no',ellip=iellip,hellip='no',interactive='yes',step=0.1,linear='no')
        else:
            iraf.ellipse(input=input_image,output=efile,dqf=mask_image,x0=xcenter_cutout,y0=ycenter_cutout,hcenter='no',recenter=0,sma0=initialr,minsma=minr,maxsma=maxr,pa=ipa,hpa='no',ellip=iellip,hellip='no',interactive='yes',step=0.1,linear='no')
    print 'Displaying isophotes from first pass.  Hit q in DS9 window to quit'
    iraf.isoexam(table=efile)
    try:
        os.system('rm tp.00*')
    except:
        print 'no tp.00? files to remove'
    iraf.isoimap(input_image,table=efile)
    # get final ellipse parameters to return to the user

    if os.path.exists('junk.txt'):
        os.remove('junk.txt')
    iraf.tprint(table=efile,pwidth='INDEF',showhdr='no', Stdout='junk.txt')
    s="awk '{print $2, $7, $9, $11, $13}' < "+"junk.txt > "+"junk2.txt"
    os.system(s)
    redo = 1
    while redo:

        infile=open('junk2.txt','r')
        print 'within redo loop'
        for line in infile:
            t=line.split()
            if float(t[0]) > myradius:
                #print 'got here',myradius, t[0]
                newellip=float(t[1])
                if newellip < .05: # min value that ellipse can handle
                    newellip=.05
                newPA=float(t[2])
                if newPA < -90:
                    newPA=newPA+180
                elif newPA > 90:
                    newPA = newPA-180
                
                newxcenter=float(t[3])
                newycenter=float(t[4])
                break
        infile.close()
        print 'ellipse values at r = ',myradius
        print '    ellip = ',newellip
        print '    PA = ',newPA
        t=raw_input('are you happy? (n = no, enter to continue) ')
        if t.find('n') >-1:
            newr=raw_input('Enter radius where ellipse best represents galaxy ')
            myradius=float(newr)        
        else:
            redo = 0
            continue

    if os.path.exists(efile_ascii):
        os.remove(efile_ascii)
    iraf.tprint(table=efile,pwidth='INDEF',showhdr='no', Stdout=efile_ascii)
    return newxcenter,newycenter,newellip,newPA

iraf.stsdas()
iraf.analysis()
iraf.isophote()
iraf.tables()
iraf.ttools()


#haimage=sys.argv[2]
#rimage=sys.argv[1]

if __name__ == '__main__':
    rimage=args.cluster+'-'+args.id+'-R.fits'
    haimage=args.cluster+'-'+args.id+'-CS.fits'

    #rimage=args.r
    #haimage=args.ha
    mask=args.cluster+'-'+args.id+'-mask.fits'
    
    # run ellipse on r image
    newxcenter,newycenter,newellip,newPA=call_measure_disk(rimage,wave_band=0,mask_image=mask,initialr=20)

    # use ellipse parameters from r-band as input for Halpha
    newxcenter,newycenter,newellip,newPA=call_measure_disk(haimage,wave_band=1,ipa=newPA,xcenter=newxcenter,ycenter=newycenter,minr=2,initialr=20,iellip = newellip,keepfixed=1,mask_image=mask)
