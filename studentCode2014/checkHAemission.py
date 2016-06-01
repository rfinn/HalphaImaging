#!/usr/bin/env python
import glob 
import ds9
import sys
import os

'''Usage - Run from command line
checkHAemission.py MKW11
'''
# open ds9  
# ask user if there is HA emission
# ask user if galaxy is a spiral or not 
# Place cutouts in designated folder 
# each galaxy will have its own folder containing Rband, Halpha with cont & w/o cont
#....................................
 
def displaycutout(prefix):
	
	d=ds9.ds9() #open new ds9 window 
	
	Rfiles = glob.glob('*R.fits') #seperate out each individual galaxy 
	#print Rfiles
	for rimage in Rfiles:
		d.set('frame delete all')
		t=rimage.split('R')
		csimage=t[0]+'Ha.fits' #Ha image w/o continuum 
		Haimage=t[0]+'Hawc.fits' #Ha image w/ continuum 
		continueflag=0
		t1=t[0].split('_') #obtain name of galaxy AGC number 
		#print rimage, os.path.isfile(rimage), os.path.isfile(csimage)
		if (os.path.isfile(rimage) & os.path.isfile(csimage) & os.path.isfile(Haimage)):
			
			s1='file new '+csimage #open Halpha cont subtracted and Rband images in ds9
			s2='file new '+rimage  #set settings for viewing in ds9
			try:
				d.set(s1)
				d.set('zscale')
				d.set(s2)
				d.set('zscale')
				
				
			except:
				print 'problem opening parent image'
                        

			while not(continueflag):
				flag=raw_input('Does this image have Halpha emission? (1 = yes, 0 = no, x = quit) \n') 
				if flag.find('x') > -1: 
					return
				try:
					flag=int(flag) # if valid command then continue to moving images
					continueflag=1
				except:
					print 'Error processing your input'
					

				spiral=raw_input('Is this a sprial galaxy? (1 = yes, 0 = no, x = quit) \n') 
				if spiral.find('x') > -1: 
					return
				try:
					spiral=int(spiral) # if valid command then continue to moving images
					continueflag=1
				except:
					print 'Error processing your input'
			
				
			if flag & spiral: #if has emission and is spiral
				soutdir=semissiondir+t1[0]+'/' 
				try:
					if os.path.isdir(soutdir):
						#print 'spiral emission directory exists'
					else:
						os.system('mkdir -p '+soutdir) #make galaxy specific folder 

					os.rename(csimage,soutdir+csimage)# move images from current folder to designated 
					os.rename(rimage,soutdir+rimage)  # output folder 
					os.rename(Haimage,soutdir+Haimage)
				except:
					print 'WARNING: problem moving image'
					
			elif flag & ~spiral: #if has emission but not spiral 
				nsoutdir=nsemissiondir+t1[0]+'/'
				try:
					if os.path.isdir(nsoutdir):
						#print 'non-spiral emission directory exists'
					else:
						os.system('mkdir -p '+nsoutdir)

					os.rename(csimage,nsoutdir+csimage)# move images 
					os.rename(rimage,nsoutdir+rimage)
					os.rename(Haimage,nsoutdir+Haimage)
				except:
					print 'WARNING: problem moving image'

			elif ~flag & spiral: #if no emission but spiral 
				
				sneoutdir=snoemissiondir+t1[0]+'/' #sne=spiral non emission (probably should change the name)
				if os.path.isdir(sneoutdir):
					#print 'non-emission spiral directory exists'
				else:
					os.system('mkdir -p '+sneoutdir)

				os.rename(csimage,sneoutdir+csimage)
				os.rename(rimage,sneoutdir+rimage)
				os.rename(Haimage,sneoutdir+Haimage)

			
			else: # if no emission or spiral 
				nsneoutdir=nsnoemissiondir+t1[0]+'/'
				if os.path.isdir(nsneoutdir):
					#print 'non-emission non spiral directory exists'
				else:
					os.system('mkdir -p '+nsneoutdir)

				os.rename(csimage,nsneoutdir+csimage)
				os.rename(rimage,nsneoutdir+rimage)
				os.rename(Haimage,nsneoutdir+Haimage)
				

	return

prefix=sys.argv[1] 

#specific output directories for the 4 differnt conditions (change to desired output location)
semissiondir='/home/share/research/LocalClusters/Halpha/cutouts/' + prefix+'/emission/spiral/'
nsemissiondir='/home/share/research/LocalClusters/Halpha/cutouts/' + prefix+'/emission/nonspiral/'
snoemissiondir = '/home/share/research/LocalClusters/Halpha/cutouts/' + prefix+'/nonemission/spiral/'
nsnoemissiondir = '/home/share/research/LocalClusters/Halpha/cutouts/' + prefix+'/nonemission/nonspiral/'

# check or create emssion/nonemssion spiral/nonspiral directories
if os.path.isdir(semissiondir): 
	print 'emission directory exists'
else:
	os.system('mkdir -p '+semissiondir)
if os.path.isdir(nsemissiondir): 
	print 'emission directory exists'
else:
	os.system('mkdir -p '+nsemissiondir)
if os.path.isdir(snoemissiondir):
	print 'non-emission directory exists'
else:
	os.system('mkdir -p '+snoemissiondir)
if os.path.isdir(nsnoemissiondir):
	print 'non-emission directory exists'
else:
	os.system('mkdir -p '+nsnoemissiondir)
	
displaycutout(prefix) #user inputed galaxy cluster name 


