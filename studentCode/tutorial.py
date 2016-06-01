# -*- coding: utf-8 -*-
"""
Spyder Editor

This temporary script file is located here:
/Users/physuser/.spyder2/.temp.py
"""
import pylab as pl
import numpy as np
import Image
def plot1():
    x=np.arange(0,2*np.pi,.1)
    y=np.sin(x)
    myerr=np.random.rand(len(y))
    pl.figure(1)
    pl.clf()
    pl.plot(x,y,'bs',label='sine')
    pl.errorbar(x,y,yerr=myerr,fmt=None)
    pl.xlabel('theta',fontsize=20)
    pl.ylabel('Amplitude',fontsize=20)
    y2=np.cos(x)
    pl.plot(x,y2,label='cosine')
    pl.axhline(y=0,color='k',ls='--')
    pl.legend()

infile1=open('dcephei.txt','r')
date=[]
mag=[]
for line in infile1:
    t=line.split()
    date.append(float(t[1]))
    mag.append(float(t[2]))
infile1.close()
date=np.array(date,'f')
mag=np.array(mag,'f')
pl.figure(2)
pl.clf()
pl.plot(date,mag,'m^')

outfile1=open('junkdata.txt','w')
for i in range(len(mag)):
    if date[i]<42000.:
        outfile1.write('%5.2f %5.2f \n'%(date[i],mag[i]))
outfile1.close()


# new plot
x=np.arange(0,2*np.pi,.1)
y1=np.sin(x)
y2=np.cos(x)
pl.figure(3)
pl.clf()
pl.subplot(2,1,1)
pl.plot(x,y1)
pl.subplot(2,1,2)
pl.plot(x,y2)
pl.savefig('myplot.png')

# new plot

iras=np.loadtxt('iras60.txt',skiprows=3)
pl.figure(4)
pl.plot(iras[:,0],iras[:,1])

# here is how you plot images
imfile='J114302.14+193859.0.cutout.jpg'
im=Image.open(imfile)
arr=np.asarray(im)
pl.imshow(arr)
pl.xticks([])
pl.yticks([])
