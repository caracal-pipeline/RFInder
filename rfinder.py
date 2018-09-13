#!/usr/bin/env python

# Import modules
import sys, string, os
import numpy as np
import yaml
import json
import glob

import logging

from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn

import warnings

#sys.path.append('/Users/maccagni/notebooks/rfinder/RFInder_modules/')
sys.path.append('/home/maccagni/programs/RFInder/RFInder_modules/')

import rfi 
import rfinder_stats as rfi_stats
import rfinder_plots as rfi_pl
import rfinder_files as rfi_files

rfi = rfi.rfi()

rfiST = rfi_stats.rfi_stats()

__author__ = "Filippo Maccagni"
__copyright__ = "Apertif"
__version__ = "1.0.0"
__email__ = "filippo.maccagni@gmail.com"
__status__ = "Development"



if not sys.warnoptions:
    warnings.simplefilter("ignore")

####################################################################################################

class rfinder:
    '''
    
    Class to investigate the RFI behaviour during observations

    '''

    C=2.99792458e5 #km/s
    HI=1.420405751e9 #Hz

    def __init__(self, file=None):
        '''
    
        Set self.logger for spectrum extraction
        Find config file
        If not specified by user load rfinder_default.yml
    
        '''

        #set self.logger
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

        self.logger.info('\t ... Reading parameter file ... \n')
        # read database here

        if file != None:
            cfg = open(file)

        else:
            #file_default = '/Users/maccagni/notebooks/RFInder/rfinder_default.yml'
            file_default = '/home/maccagni/programs/RFInder/rfinder_default.yml'

            cfg = open(file_default) 

        self.cfg_par = yaml.load(cfg)
    
        self.set_dirs()


    def enable_task(self,config,task):

        a = config.get(task, False)
        if a:
            return a['enable']
        else:
            False

    def set_dirs(self):
        '''
     
        Sets directory strucure and filenames
        Creates directory abs/ in basedir+beam and subdirectories spec/ and plot/

        OUTPUT:
            tab : table of catalog
            flux : flux of continuum sources abouve set threshold
     
        '''

        key = 'general'

        self.workdir  = self.cfg_par[key].get('workdir', None)
        self.msfile = self.workdir + self.cfg_par[key].get('msname', None)[0]
        self.cfg_par[key]['msfullpath'] = self.msfile      

        self.rfidir  = self.workdir+'rfi/'
        self.cfg_par[key]['rfidir'] = self.rfidir

        self.rfifile = self.rfidir+'rfi_flagged_vis.MS'
        self.rfi_freq_base = self.rfidir+'freq_base.fits'
        self.rfimsfile = self.rfidir+'rfi_flagged.MS'
        self.tabledir = self.rfidir+'tables/'
        self.cfg_par[key]['tabledir'] = self.tabledir
        self.rfi_table = self.tabledir+'rfi_table.fits'
    
        self.rfiplotdir = self.rfidir+'plots/'
        self.cfg_par[key]['plotdir'] = self.rfiplotdir 

 


        if os.path.exists(self.rfidir) == False:
             os.makedirs(self.rfidir)           

        if os.path.exists(self.tabledir) == False:
             os.makedirs(self.tabledir)

        if os.path.exists(self.rfiplotdir) == False:
             os.makedirs(self.rfiplotdir)

        if self.cfg_par['rfi']['chunks']['time_enable'] == True:

            self.rfitimedir = self.rfidir+'time_chunks/'
            self.cfg_par[key]['rfitimedir'] = self.rfitimedir

            if os.path.exists(self.rfitimedir) == False:
                 os.makedirs(self.rfitimedir)

            self.timetabledir = self.tabledir+'time_chunks/'
            self.cfg_par[key]['timetabledir'] = self.timetabledir

            if os.path.exists(self.timetabledir) == False:
                 os.makedirs(self.timetabledir)

            self.timeplotdir = self.rfiplotdir+'time_chunks/'
            self.cfg_par[key]['timeplotdir'] = self.timeplotdir

            if os.path.exists(self.timeplotdir) == False:
                 os.makedirs(self.timeplotdir)

            self.altazplotdir = self.rfidir+'plots/altaz/'
            self.cfg_par[key]['altazplotdir'] = self.altazplotdir        

            if os.path.exists(self.altazplotdir) == False:
                 os.makedirs(self.altazplotdir)

    def go(self,cfg_par):
        '''
        Automated pipeline to extract spectra from each continuum source in a given field.
        If cfg_par['rfi'] is enabled 
            Executes the whole spectrum extraction process as follows:
            1: load_from_ms
            2: baselines_from_ms 
            3: priors_flag
            4: rfi_flag
        If cfg_par['plots'] is enabled
            1: 2d plot of RFI flagged by frequency and baseline lenght (plot_rfi_im)
            2: 1d plot of RFI flagged by frequency channel (baselines_from_ms)
            3: 1d plot of noise increase by frequency channel (for long and short baselines) (priors_flag)
        If cfg_par['beam_shape'] is enabled
            1: create FLAG column in MS file (rfi_flag)
            2: determine psf using wsclean (make_psf)
        '''

        # cont_sources
        task = 'rfi'
        self.logger.info(" ------ STARTING RFI analysis ------\n")

        if self.enable_task(self.cfg_par,task)==True:
            
            if self.cfg_par[task]['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.info(" ------ Working on time chunks ------\n")

                for i in xrange(0,len(times)-1):
                    timez = [times[i],times[i+1]] 
                    self.logger.info((" ------ Working on chunk #{0:d}: {1:d} minutes after begin of observation ------\n").format(i,int(cfg_par['rfi']['chunks']['time_step']*i)))

                    rfi.load_from_ms(self.cfg_par,timez)
                    #sort visibilities by baseline lenght
                    rfi.baselines_from_ms(self.cfg_par)

                    #flag bad antennas (from configuration file)
                    datas = rfi.priors_flag(self.cfg_par)

                    #find rfi above threshold
                    rfi.find_rfi(datas,self.cfg_par,i)

                    rfi_files.rfi_frequency(self.cfg_par,i)
                    rfi_pl.plot_rfi_imshow(self.cfg_par,i)
                    rfi_pl.plot_noise_frequency(self.cfg_par,i)

                rfi_pl.plot_altaz(self.cfg_par,i)

                self.logger.info(" ------ End of RFI analysis on time chunks ------\n")

            else:

                rfi.load_from_ms(self.cfg_par,0)
                #determine alt/az

                self.logger.info("---- MSfile Loaded -----\n")
                rfi.baselines_from_ms(self.cfg_par)
                self.logger.info("---- Dataset sorted by baseline lenght ----\n")
                datas = rfi.priors_flag(self.cfg_par)
                self.logger.info("---- Bad antennas and autocorrelations flagged ----\n")
                rfi.find_rfi(datas,self.cfg_par,-1)
        
        task = 'plots'
        if self.enable_task(self.cfg_par,task) == True:
            
            if self.cfg_par['rfi']['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.info(" ------ Working on time chunks ------\n")

                for i in xrange(0,len(times)-1):
                    timez = [times[i],times[i+1]]            

                    if self.enable_task(self.cfg_par,'rfi')==False:

                        rfi.load_from_ms(self.cfg_par,timez)
                        self.logger.info("---- MSfile Loaded -----\n")                
                        rfi.baselines_from_ms(self.cfg_par)
                        self.logger.info("---- Dataset sorted by baseline lenght ----\n")
            
                    rfi_pl.plot_rfi_imshow(self.cfg_par,i)
                    self.logger.info("---- RFI in 2D plotted ----\n")
                    rfi_files.rfi_frequency(self.cfg_par,i)
                    self.logger.info("---- RFI saved to table ----\n")
                    rfi_pl.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['long_short'] = True
                    self.cfg_par['plots']['plot_noise'] = 'noise_factor'
                    rfi_pl.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_noise'] = 'noise'
                    rfi_pl.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['long_short'] = False  
                    rfi_pl.plot_noise_frequency(self.cfg_par,i)            
                    self.logger.info("---- RFI in 1D plotted ----\n")
                rfi_pl.plot_altaz(self.cfg_par,i)
                self.logger.info("---- RFI in ALT/AZ plotted ----\n")

            else:

                rfi_pl.plot_rfi_imshow(self.cfg_par,-1)
                self.logger.info("---- RFI in 2D plotted ----\n")
                rfi_files.rfi_frequency(self.cfg_par,-1)
                self.logger.info("---- RFI saved to table ----\n")
                rfi_pl.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['long_short'] = True
                self.cfg_par['plots']['plot_noise'] = 'noise_factor'
                rfi_pl.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_noise'] = 'noise'
                rfi_pl.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['long_short'] = False  
                rfi_pl.plot_noise_frequency(self.cfg_par,-1)            
                self.logger.info("---- RFI in 1D plotted ----\n")

        task = 'beam_shape'
        if self.enable_task(self.cfg_par,task) == True:
            rfi.rfi_flag(self.cfg_par)
            rfi_beam.make_psf(self.cfg_par)
        self.logger.info(" ------ End of RFI analysis ------ \n\n")

    
        return 0






