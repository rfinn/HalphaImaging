

class cluster:

    
    def just_get_images24(self,i, keepimages=1,make_mask_flag=0,review_mask_flag=0):
        # getimages:
        # set = 0 to not get images (if they already exist)
        # set = 1 if sdss images are needed
        #
        # this subroutine gets a cutouot from the mips mosaic,
        #    parses relevant input for galfit,
        #    and runs galfit three times.
        #        #
        # keepimages = 1 to keep MIPS cutout images
        # keepimages = 0 to delete MIPS cutout images (better if disk space is tight)
        quitflag=0
        ##############################################
        # MAG ZP
        ##############################################
        # define mag zp for 24um image scans
        flux_zp_AB = 3631. # in Jy
        flux_zp_Vega = 7.17 # in Jy
        flux_zp=flux_zp_AB
        
        # conversion from image units of MJ/sr to micro-Jy (1 sq arcsec = 2.3504e-11 sr)
        conv_MJysr_uJy = 23.5045*(2.45**2)
        magzp=2.5*log10(flux_zp*1.e6/conv_MJysr_uJy)

        ##############################################
        # SOME SETUP COMMANDS
        ##############################################
        # open ds9 display
        d=ds9.ds9()
        working_dir=homedir+'research/LocalClusters/GalfitAnalysis/'+self.prefix+'/24um/'
        s='mkdir -p '+working_dir
        os.system(s)
        os.chdir(working_dir)
        d.set('cd '+working_dir)
        # clean up any previously created galfit files
        os.system('rm fit.log')
        os.system('rm galfit.??')

        magzp=magzp+5.


        ra=self.ra[i]
        dec=self.dec[i]
        s='echo %f %f > '%(ra,dec)
        t=s+working_dir+'incoords'
        os.system(t)
        print 'deleting outcoords if it exists'
        output_coords=working_dir+'outcoords'
        if os.path.exists(output_coords):
            os.remove(output_coords)
        input_coords=working_dir+'incoords'
        print input_coords,os.getcwd()
        iraf.imcoords.wcsctran(image=self.mosaic24,input=input_coords,output=output_coords,inwcs='world',outwcs='logical',verbose='no')
        coordout=atpy.Table(output_coords,type='ascii')
        col=coordout['col1']
        row=coordout['col2']
        # convert arrays to numbers
        col=col[0]
        row=row[0]


        ##############################################
        # GET CUTOUT OF GALAXY AND UNCERTAINTY IMAGE
        ##############################################
        radius=self.n.PETROTH90[i]
        #radius=2.*self.n.SERSIC_TH50[i]
        if radius > min_cutout_radius:
            PETROTH90_pixels =radius/mipspixelscale
            if radius> (max_cutout_radius):
                    PETROTH90_pixels =max_cutout_radius/mipspixelscale

        else:
            PETROTH90_pixels =(min_cutout_radius)/mipspixelscale

        # make a cutout of the galaxy (easier than trying to run sextractor on full mosaic
        # need x and y pixel values of the galaxy
        xmin=col-mult_petro90*PETROTH90_pixels
        xmax=col+mult_petro90*PETROTH90_pixels
        ymin=row-mult_petro90*PETROTH90_pixels
        ymax=row+mult_petro90*PETROTH90_pixels

        # get image dimensions of 24um mosaic

        iraf.imgets(image=self.mosaic24,param='naxis1')#get RA of image
	image_xmax=float(iraf.imgets.value)
        iraf.imgets(image=self.mosaic24,param='naxis2')#get RA of image
	image_ymax=float(iraf.imgets.value)
        
        # check that cutout region is not outside bounds of image

        if (xmax > image_xmax):
            print 'WARNING: fit region extends beyond image boundary (image size = %5.2f, xmaxfit = %5.2f)'%(image_xmax,xmax)
            print self.prefix, self.n.NSAID[i], ': setting xmax to max x of image'
            xmax=image_xmax

        if (ymax > image_ymax):
            print 'WARNING: fit region extends beyond image boundary (image size = %5.2f, xmaxfit = %5.2f)'%(image_ymax,ymax)
            print self.prefix, self.n.NSAID[i], ': setting ymax to max y of image'
            ymax=image_ymax

        if (xmin < 1):
            print 'WARNING: fit region extends beyond image boundary (image size = %5.2f, xminfit = %5.2f)'%(image_xmax,xmin)
            print self.prefix, self.n.NSAID[i], ': setting xmin to 1'
            xmin=1

        if (ymin < 1):
            print 'WARNING: fit region extends beyond image boundary (image size = %5.2f, yminfit = %5.2f)'%(image_xmax,ymin)
            print self.prefix, self.n.NSAID[i], ': setting ymin to 1'
            ymin=1
        
        sex_image=self.prefix+'-'+str(self.n.NSAID[i])+'-'+'galfit-cutout24.fits'
        unc_image=self.prefix+'-'+str(self.n.NSAID[i])+'-'+'galfit-cutout-unc24.fits'

        try:
            s=self.mosaic24+'[%i:%i,%i:%i] '%(int(round(xmin)),int(round(xmax)),int(round(ymin)),int(round(ymax)))
            iraf.imcopy(s,working_dir+sex_image)
            s=self.mosaic24unc+'[%i:%i,%i:%i] '%(int(round(xmin)),int(round(xmax)),int(round(ymin)),int(round(ymax)))
            iraf.imcopy(s,working_dir+unc_image)

        except:
            print 'Warning:  Problem creating cutout image, probably b/c it already exists'
            print '   going to delete existing cutout image and try imcopy again'
            iraf.imdel(sex_image)
            iraf.imdel(unc_image)
            s=self.mosaic24+'[%i:%i,%i:%i] '%(int(round(xmin)),int(round(xmax)),int(round(ymin)),int(round(ymax)))
            iraf.imcopy(s,working_dir+sex_image)
            s=self.mosaic24unc+'[%i:%i,%i:%i] '%(int(round(xmin)),int(round(xmax)),int(round(ymin)),int(round(ymax)))
            iraf.imcopy(s,working_dir+unc_image)

        ##############################################
        # RUN SEXTRACTOR AND MAKE OBJECT MASK
        ##############################################

        # run sextractor to generate a list of objects in the image and generate 'segmentation image'
        os.system('cp ~/research/LocalClusters/sextractor/default.sex.24um.galfit .')
        os.system('cp ~/research/LocalClusters/sextractor/default.param .')
        os.system('cp ~/research/LocalClusters/sextractor/default.conv .')
        os.system('cp ~/research/LocalClusters/sextractor/default.nnw .')
        os.system('sex %s -c default.sex.24um.galfit'%(sex_image))
        # convert segmentation image to object mask by replacing the object ID of target galaxy with zeros
        #   parse sextractor output to get x,y coords of objects        
        sexout=atpy.Table('test.cat',type='ascii')
        sexnumber=sexout['col1']
        xsex=sexout['col2']
        ysex=sexout['col3']
        xcenter=1.*(xmax-xmin)/2.
        ycenter=1.*(ymax-ymin)/2.
        dist=sqrt((ycenter-ysex)**2+(xcenter-xsex)**2)

        #   find object ID
        objIndex=where(dist == min(dist))
        objNumber=sexnumber[objIndex]
        objNumber=objNumber[0] # not sure why above line returns a list
        print 'object number = ',objNumber
        if make_mask_flag:
            print working_dir+self.prefix+'-'+str(self.n.NSAID[i])+'-'+'galfit-mask24.fits'
            mask_image=working_dir+self.prefix+'-'+str(self.n.NSAID[i])+'-'+'galfit-mask24.fits'

            print mask_image
            if os.path.exists(mask_image):
                iraf.imdel(mask_image)
            try:
                iraf.imcopy('segmentation.fits',mask_image)
            except:
                iraf.imcopy('segmentation.fits[1]',mask_image)

        #   use iraf imexam to replace object ID values with zero
            iraf.imreplace(mask_image,value=0,lower=objNumber-.5,upper=objNumber+.5)

        while review_mask_flag:
            d.set('frame delete all')
            s='file new '+homedir+'research/LocalClusters/GalfitAnalysis/'+self.prefix+'/NSA/'+self.prefix+'-'+str(self.n.NSAID[i])+'-1Comp-galfit-out.fits[1]'
            try:
                d.set(s)
                d.set('zoom to fit')
            except:
                print "couldn't access: ",s
            s='file new '+sex_image
            d.set(s)
            d.set('zoom to fit')
            s='file new '+mask_image
            d.set(s)
            d.set('zoom to fit')
            flag=raw_input('edit the mask? \n')
            if flag.find('n') > -1:
                review_mask_flag = 0
            elif flag.find('y') > -1:
                editflag=int(raw_input('enter 0 to subtract an object, 1 to add an object'))
                if editflag == 0:
                    objid=float(raw_input('enter object id to remove from mask'))
                    iraf.imreplace(mask_image,value=0,lower=objid-.5,upper=objid+.5)
                elif editflag == 1:
                    print 'entering imedit'
                    a,b=runimedit(mask_image)
            elif flag.find('q') > -1:
                quitflag=1
                return quitflag

        # get best x,y positions from sextractor
        #    - subtract 1 b/c sextractor catalog starts at one, but python is zero-indexed
        xcenter=xsex[objNumber-1]
        ycenter=ysex[objNumber-1]

        self.xcenter_cutout24[i]=(xcenter)
        self.ycenter_cutout24[i]=(ycenter)
        return quitflag


