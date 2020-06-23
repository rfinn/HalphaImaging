'''
Python version of Becky Koopmann's idl program fluxconc_groups_magcut_notype_ex.pro


python version by Rose Finn
written Oct 22, 2019

'''

#!/usr/bin/env python
from astropy.io import ascii, fits
from astropy.table import Table
import os
import numpy as np
from matplotlib import pyplot as plt
import glob

######################################################
# change this to directory where date files are
######################################################

data_dir = os.getenv('HOME')+'/research/HalphaGroups/Becky-sSFR-concentration/'
class table_tools():
    def read_tables(self):
        # read in virgo r
        #self.rtab = ascii.read(rfile,delimiter='\s',names = colnames_allr, guess=False, fast_reader=False, data_start=0)
        mycols = np.arange(len(self.colnames_r)).tolist()
        rtab = np.loadtxt(self.rfile,unpack=True,usecols = mycols)
        self.rtab = Table(np.array(rtab).T, names = self.colnames_r)
        #self.rtab = ascii.read(rfile,delimiter='\s',names = colnames_allr, guess=False, fast_reader=False, data_start=0)
        
        # read in virgo ha
        mycols = np.arange(len(self.colnames_h)).tolist()
        htab = np.loadtxt(self.hfile,unpack=True,usecols = mycols)
        self.htab = Table(np.array(htab).T, names = self.colnames_h)

        #self.htab = ascii.read(hfile,delimiter='\s',names = colnames_allh, guess=False, fast_reader=False,data_start=0)

class get_sfr():
    def calc_sfr(self):
        self.haflux=np.log10(self.htab['tefl']*1e-18)
        self.rflux=(-(self.rtab['tmag24']+13.945)/2.5)
        self.nhaflux=self.haflux-self.rflux

        self.haflux30=np.log10((self.htab['inn30fl']*1e-18)/0.9)
        self.rflux30=np.log10(self.rtab['c30']*10**(-(self.rtab['tmag24']+13.945)/2.5))
        self.nhaflux30=self.haflux30-self.rflux30
        
        self.haflux70 = np.log10(10.**self.haflux - 10.**self.haflux30)
        self.rflux70 = np.log10(10.**self.rflux - 10.**self.rflux30)
        self.nhaflux70 = self.haflux70 - self.rflux70
        
        self.set_low_flag()
    def set_low_flag(self):
        self.nhafluxlow=-2.935*self.rtab['c30']-1.013
        self.low_flag = self.nhaflux < self.nhafluxlow
    def calc_sfr_group(self):
        self.haflux=np.log10((self.htab['tefl']*1e-18)/0.9)
        self.rflux=(-(self.rtab['tmag24']+13.945)/2.5)
        self.nhaf=self.haflux-self.rflux
        self.nhaflow=-2.935*self.rtab['c30']-1.013


    def plot_ssfr_conc(self,ax1=None,ax2=None,ax3=None,plotsingle=False,filament=False):
        if plotsingle:
            fig = plt.figure()
            ax1 = fig.subplot(1,3,1)
            ax3 = fig.subplot(1,3,2)
            ax3 = fig.subplot(1,3,3)


        all_ax = [ax1, ax2, ax3]
        all_yvar = [self.nhaflux, self.nhaflux30, self.nhaflux70]

        for i,ax in enumerate(all_ax):
            print(i)
            self.add_lines(ax,i)
            ax.plot(self.rtab['c30'],all_yvar[i],'ko')
            #if not filament:
            ax.plot(self.rtab['c30'][self.low_flag],all_yvar[i][self.low_flag],'bo')
            self.set_limits(ax)
            if i >-1:
                print('printing name')
                ax.text(.2,-4,self.name,horizontalalignment='left',fontsize=12)
        
    def add_lines(self,ax,i):

        if i == 0:
            ax.plot([0.2,0.66],[-1.6,-2.95],'c--')
            ax.plot([0.2,0.66],[-1.4,-1.7],'k--')
        if i == 2:
            ax.plot([0.2,0.64],[-1.5,-3.1],'k--')
            ax.plot([0.2,0.6],[-1.3,-1.5],'k--')
        if i == 1:
            ax.plot([0.2,0.64],[-2.1,-3.4],'k--')
            ax.plot([0.2,0.6],[-1.7,-1.7],'k--')
        #ax.plot([0.2,0.66],[-1.6,-2.95],'c--')
    def set_limits(self,ax):
        ax.axis([.18,.72,-4.13,-.99])
        
class filament(get_sfr):
    def __init__(self,pfile):
        self.name = 'Filament'
        self.rtab = fits.getdata(pfile)
        self.nhaflux = np.log10(self.rtab.HF_R24) - np.log10(self.rtab.F_R24)
        self.nhaflux30 = np.log10(self.rtab.HF_30R24) - np.log10(self.rtab.F_30R24)
        self.nhaflux70 = np.log10(self.rtab.HF_R24 - self.rtab.HF_30R24) - np.log10(self.rtab.F_R24 - self.rtab.F_30R24)
        self.set_low_flag()
        #self.plot_ssfr_conc()
               
class galaxy_group(get_sfr):
    def __init__(self,rfile, hfile):
        self.name='Groups'
        # read input file
        colnames_r = ['gal','inc','dist','r24','r24e',\
                      'r25','r25e','r26','r26e','lmag',\
                      'tmag24','tmag25', 'c30','c7525','reff24']
        colnames_h = ['galha','inc','dist','rl','r95',\
                      'rr9524','tefl','r17','rr1724','r17fl',\
                      'cnfl','inn30fl','totr24flux','c30ha','c7525ha',\
                      'reffhar24']
        self.rtab = ascii.read(rfile,delimiter='\s',names = colnames_r, guess=False, fast_reader=False)
        self.htab = ascii.read(hfile,delimiter='\s',names = colnames_h, guess=False, fast_reader=False)
        self.calc_sfr()
        
class virgo_cluster(get_sfr,table_tools):

    '''
    readcol,'virgo_r.sdat',vgal,vtype,r24,r24e,r25,r25e,lmag, $
    vtmag24,tmag25, vc30,c7525,reff24,ir
    
    readcol,'virgo_ha.sdat',galha,rl,r95,rr9524,vtefl, $
    r17,rr1724,r17fl,cnfl,vinn30fl,c30ha,c30har30,c7525ha,reffha,reffhar24,iha
    '''
    def __init__(self,rfile,hfile):
        self.name='Virgo'
        self.rfile = rfile
        self.hfile = hfile
        self.colnames_r = ['gal','type','r24','r24e',\
                      'r25','r25e','lmag',\
                      'tmag24','tmag25', 'c30','c7525','reff24','ir']
        self.colnames_h = ['gal','rl','r95',\
                      'rr9524','tefl','r17','rr1724','r17fl',\
                      'cnfl','inn30fl','c30ha','c30har30','c7525ha',\
                      'reffha','reffhar24','iha']




        self.read_tables()
        # calculate sSFR
        self.calc_sfr()

        pass
class field(get_sfr,table_tools):
    '''
    readcol,'iso_r.sdat',igal,itype,r24,r24e,r25,r25e,lmag, $
    itmag24,tmag25, ic30,c7525,reff24,ir

    readcol,'iso_ha.sdat',galha,rl,r95,rr9524,itefl, $
    r17,rr1724,r17fl,cnfl,iinn30fl,c30ha,c30har30,c7525ha,reffha,reffhar24,iha
    '''
    def __init__(self,rfile, hfile):
        self.name='Field'
        self.rfile = rfile
        self.hfile = hfile
        self.colnames_r = ['gal','type','r24','r24e',\
                      'r25','r25e','lmag',\
                      'tmag24','tmag25', 'c30','c7525','reff24','ir']
        self.colnames_h = ['gal','rl','r95',\
                      'rr9524','tefl','r17','rr1724','r17fl',\
                      'cnfl','inn30fl','c30ha','c30har30','c7525ha',\
                      'reffha','reffhar24','iha']
        # read in virgo r
        #self.rtab = ascii.read(rfile,delimiter='\s',names = colnames_r, guess=False, fast_reader=False)
        # read in virgo ha
        #self.htab = ascii.read(hfile,delimiter='\s',names = colnames_h, guess=False, fast_reader=False)

        self.read_tables()
        # calculate sSFR
        self.calc_sfr()

        pass




rfile = data_dir+'virgo_r_clean.sdat'
hfile = data_dir+'virgo_ha.sdat'
v = virgo_cluster(rfile,hfile)

# isolated galaxies
rfile = data_dir+'iso_r_clean.sdat'
hfile = data_dir+'iso_ha_clean.sdat'
f = field(rfile,hfile)


def plotit():
    nrow = 3
    ncol = 3
    fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_subplot(nrow,ncol,1)
    plt.title('GLOBAL',fontsize=16)
    ax2 = fig.add_subplot(nrow,ncol,2)
    plt.title('INNER 30%',fontsize=16)
    ax3 = fig.add_subplot(nrow,ncol,3)
    plt.title('OUTER 70%',fontsize=16)
    f.plot_ssfr_conc(ax1=ax1,ax2=ax2,ax3=ax3)
    ax1 = fig.add_subplot(nrow,ncol,4)
    plt.ylabel(r'$\log_{10}(F_{H \alpha}/F_{R})$',fontsize=20)
    ax2 = fig.add_subplot(nrow,ncol,5)

    ax3 = fig.add_subplot(nrow,ncol,6)
    v.plot_ssfr_conc(ax1=ax1,ax2=ax2,ax3=ax3)
    ax1 = fig.add_subplot(nrow,ncol,7)
    ax2 = fig.add_subplot(nrow,ncol,8)
    plt.xlabel('$C30$',fontsize=20)
    ax3 = fig.add_subplot(nrow,ncol,9)

    ggroups = ['NRGb027','NRGb032','NRGb151','NRGb282','NRGb317','NRGb331','NRGs272']

    for g in ggroups:
        try:
            print(g)
            rfile = data_dir+g+'r.sdat'
            hfile = data_dir+g+'ha.sdat'
            gg = galaxy_group(rfile,hfile)

            gg.plot_ssfr_conc(ax1=ax1,ax2=ax2,ax3=ax3)
        except:
            print('problem with group ',g)


    plt.savefig('ssfr-c30.png')


def plotfilaments():
    nrow = 4
    ncol = 3
    fig = plt.figure(figsize=(10,8))
    ax1 = fig.add_subplot(nrow,ncol,1)
    plt.title('GLOBAL',fontsize=16)
    ax2 = fig.add_subplot(nrow,ncol,2)
    plt.title('INNER 30%',fontsize=16)
    ax3 = fig.add_subplot(nrow,ncol,3)
    plt.title('OUTER 70%',fontsize=16)
    f.plot_ssfr_conc(ax1=ax1,ax2=ax2,ax3=ax3)
    ax1 = fig.add_subplot(nrow,ncol,4)
    plt.ylabel(r'$\log_{10}(F_{H \alpha}/F_{R})$',fontsize=20)
    ax2 = fig.add_subplot(nrow,ncol,5)

    ax3 = fig.add_subplot(nrow,ncol,6)
    v.plot_ssfr_conc(ax1=ax1,ax2=ax2,ax3=ax3)
    
    ax1 = fig.add_subplot(nrow,ncol,7)
    ax2 = fig.add_subplot(nrow,ncol,8)
    #plt.xlabel('$C30$',fontsize=20)
    ax3 = fig.add_subplot(nrow,ncol,9)

    ggroups = ['NRGb027','NRGb032','NRGb151','NRGb282','NRGb317','NRGb331','NRGs272']

    for g in ggroups:
        try:
            print(g)
            rfile = data_dir+g+'r.sdat'
            hfile = data_dir+g+'ha.sdat'
            gg = galaxy_group(rfile,hfile)

            gg.plot_ssfr_conc(ax1=ax1,ax2=ax2,ax3=ax3)
        except:
            print('problem with group ',g)

    ax1 = fig.add_subplot(nrow,ncol,10)
    ax2 = fig.add_subplot(nrow,ncol,11)
    plt.xlabel('$C30$',fontsize=20)
    ax3 = fig.add_subplot(nrow,ncol,12)


    ffiles = glob.glob('/Users/rfinn/research/VirgoFilaments/Halpha/results/2017/v17p*rfinn*.fits')
    for p in ffiles:
        print(p)
        fil = filament(p)
        fil.plot_ssfr_conc(ax1=ax1,ax2=ax2,ax3=ax3,filament=True)
        #except:
        #    print('problem with filament ',p) 
    plt.savefig('ssfr-c30-filaments.png')
