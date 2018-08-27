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

sys.path.append('/Users/maccagni/notebooks/rfinder/RFInder_modules/')
import rfi 
import rfinder_beam as rfi_beam
import rfinder_plots as rfi_plots


rfi = rfi.rfi()

__author__ = "Filippo Maccagni"
__copyright__ = "RadioLife"
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
    
        Set logger for spectrum extraction
        Find config file
        If not specified by user load rfinder_default.yml
    
        '''

        #set logger
        self.logger = logging.getLogger('RFI_general')


        file_default = '/Users/maccagni/notebooks/rfinder/rfinder_default.yml'

        if file != None:
            cfg = open(file)

        else:
            cfg = open(file_default) 

        self.cfg_par = yaml.load(cfg)
    
        self.set_dirs()

        return

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
        self.msfile = self.workdir + self.cfg_par[key].get('msname', None)
        self.cfg_par[key]['msname'] = self.msname       

        self.rfidir  = self.workdir+'rfi/'
        self.cfg_par[key]['rfidir'] = self.rfidir

        self.rfifile = self.rfidir+'rfi_flagged_vis.MS'
        self.rfi_freq_base = self.rfidir+'freq_base.fits'
        self.rfimsfile = self.rfidir+'rfi_flagged.MS'
        self.rfi_table = self.rfidir+'rfi_table.fits'
    
        self.rfiplotdir = self.rfidir+'plot/'
        self.cfg_par[key]['plotdir'] = self.plotdir 

        self.rfi_freq_plot = self.rfiplotdir+'freq_fri.png'
        self.rfi_freq_base_plot = self.rfiplotdir+'freq_base.png'



        if os.path.exists(self.rfidir) == False:
             os.makedirs(self.rfidir)           

        if os.path.exists(self.plotdir) == False:
             os.makedirs(self.plotdir)



    def go(self,cfg_par):
        '''
        Automated pipeline to extract spectra from each continuum source in a given field.

        Executes the whole spectrum extraction process as follows:
        1: load_from_ms
        2: baselines_from_ms 
        3: priors_flag
        4: rfi_flag
        '''
        catalog_table = str(self.cfg_par['general'].get('absdir')) + 'cat_src_sharpener.txt'

        # cont_sources
        task = 'rfi'
        self.logger.info("########## STARTING RFI analysis ##########")

        if self.enable_task(self.cfg_par,task)==True:
            rfi.load_from_ms(self.cfg_par,catalog_table)
            rfi.baselines_from_ms(self.cfg_par)
            self.logger.info("---- Baselines written -----")
            rfi.priors_flag(self.cfg_par)
            self.logger.info("---- Bad antennas and autocorrelations flagged ----")
            rfi.rfi_flag(self.cfg_par)
        self.logger.info("---- RFI found ----")
        
        task = 'plots'
        if self.enable_task(self.cfg_par,task) == True:
            rfi_plots.plot_rfi_im(self.cfg_par)

        task = 'beam_shape'
        if self.enable_task(self.cfg_par,task) == True:
            rfi_beam.make_psf(self.cfg_par)

    
        return 0






