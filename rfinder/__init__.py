# Import modules
import os
import sys
import string
import numpy as np
import yaml
import json
import glob

import logging

from astropy.io import fits, ascii
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.table import Table, Column, MaskedColumn

import warnings

# get rfinder install directory
RFINDER_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(os.path.join(RFINDER_PATH, 'rfinder'))

import rfi
import rfinder_stats as rfi_stats
import rfinder_plots as rfi_plots
import rfinder_files as rfi_files

rfi = rfi.rfi()
rfiST = rfi_stats.rfi_stats()
rfiPL = rfi_plots.rfi_plots()

DEFAULT_CONFIG = 'rfinder_default.yml'

if not sys.warnoptions:
    warnings.simplefilter("ignore")

####################################################################################################


class rfinder:
    '''

    Class to investigate the RFI behaviour during observations

    '''

    C = 2.99792458e5  # km/s
    HI = 1.420405751e9  # Hz

    def __init__(self, file=None):
        '''

        Set self.logger for spectrum extraction
        Find config file
        If not specified by user load rfinder_default.yml

        '''

        # set self.logger
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

        self.logger.info('\t ... Reading parameter file ... \n')
        # read database here

        if file is not None:
            cfg = open(file)

        else:
            file_default = os.path.join(RFINDER_PATH, DEFAULT_CONFIG)

            cfg = open(file_default)

        self.cfg_par = yaml.load(cfg)
        self.cfg_par['general']['template_folder'] = os.path.join(RFINDER_PATH,'rfinder/templates')
        self.set_dirs()

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

        self.outdir  = self.cfg_par[key].get('outdir', None)
        self.rfidir  = self.outdir+'rfi_'+self.cfg_par['rfi']['polarization']+'/'
        self.cfg_par[key]['rfidir'] = self.rfidir
        self.rfifile = self.rfidir+'rfi_flagged_vis.MS'
        self.rfi_freq_base = self.rfidir+'freq_base.fits'
        self.rfimsfile = self.rfidir+'rfi_flagged.MS'
        self.tabledir = self.rfidir+'tables/'
        self.cfg_par[key]['tabledir'] = self.tabledir
        self.rfi_table = self.tabledir+'rfi_table.fits'
    
        self.rfiplotdir = self.rfidir+'plots/'
        self.cfg_par[key]['plotdir'] = self.rfiplotdir 

 
        self.moviedir = self.rfidir+'plots/movies/'
        self.cfg_par[key]['moviedir'] = self.moviedir        

        if os.path.exists(self.moviedir) == False:
             os.makedirs(self.moviedir)

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

            timeplotdir_tmp = self.rfiplotdir+'time_chunks/'
            
            if os.path.exists(timeplotdir_tmp) == False:
                 os.makedirs(timeplotdir_tmp)

            self.timeplotdir1d = timeplotdir_tmp+'1D/'
            self.cfg_par[key]['timeplotdir1D'] = self.timeplotdir1d

            if os.path.exists(self.timeplotdir1d) == False:
                 os.makedirs(self.timeplotdir1d)

            self.timeplotdir2d = timeplotdir_tmp+'2D/'
            self.cfg_par[key]['timeplotdir2D'] = self.timeplotdir2d

            if os.path.exists(self.timeplotdir2d) == False:
                 os.makedirs(self.timeplotdir2d)

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
        If cfg_par['rfi']['chunks']['time_chunks'] is enabled
            1: executes 'rfi' and 'plots' procedure dividing the observation in time-steps given by cfg_par['rfi']['chunks']['time_step']
            2: collects the info about the % of RFI for each time step in Alt/Az plots 
        If cfg_par['beam_shape'] is enabled
            1: create FLAG column in MS file (rfi_flag)
            2: determine psf using wsclean (make_psf)
        '''

        # cont_sources
        task = 'rfi'
        self.logger.info(" ------ STARTING RFI analysis ------\n")

        if self.cfg_par[task]['rfi_enable']==True:
            
            if self.cfg_par[task]['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.info(" ------ Working on time chunks ------\n")

                for i in xrange(0,len(times)-1):
                    timez = [times[i],times[i+1]] 
                    
                    #time chunk properties
                    time_delta = float(self.cfg_par['rfi']['chunks']['time_step'])*i
                    time_del = TimeDelta(time_delta*60., format='sec')
                    time_delta_plus = TimeDelta(float(self.cfg_par['rfi']['chunks']['time_step'])*60., format='sec')
                    start = self.cfg_par['rfi']['startdate']+time_del
                    end = start+time_delta_plus

                    self.logger.info((" ------ Working on chunk #{0:d}:").format(i))
                    self.logger.info(("\t \t between {0:%d}{0:%b}{0:%y}: {0:%H}:{0:%M} - {1:%H}:{1:%M}").format(start.datetime,end.datetime))

                    result = rfi.load_from_ms(self.cfg_par,timez)
                    self.logger.info("---- MSfile Loaded -----\n")

                    #sort visibilities by baseline lenght
                    if result != 1:
                        rfi.baselines_from_ms(self.cfg_par)
                        self.logger.info("---- Dataset sorted by baseline lenght ----\n")

                        #flag bad antennas (from configuration file)
                        datas = rfi.priors_flag(self.cfg_par)
                        self.logger.info("---- Bad antennas and autocorrelations flagged ----\n")

                        #find rfi above threshold
                        rfi.find_rfi(datas,self.cfg_par,i)
                        self.logger.info(" ------  RFI found  ------\n")

                        rfi_files.rfi_frequency(self.cfg_par,i)
                        self.logger.info("---- RFI saved to table ----\n")
                    else:
                        self.logger.warning("---- This chunk is empty ----\n")
                        continue

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
                self.logger.info(" ------  RFI found  ------\n")
                rfi_files.rfi_frequency(self.cfg_par,-1)
                self.logger.info("---- RFI saved to table ----\n")
                self.logger.info(" ------ End of RFI analysis  ------\n")
      
        task = 'plots'
        if self.cfg_par[task]['plot_enable']==True:
            
            if self.cfg_par['rfi']['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.info(" ------ Plotting on time chunks ------\n")

                for i in xrange(0,len(times)-1):
#                for i in xrange(0,2):

                    timez = [times[i],times[i+1]]            
                    
                    #time chunk properties
                    time_delta = float(self.cfg_par['rfi']['chunks']['time_step'])*i
                    time_del = TimeDelta(time_delta*60., format='sec')
                    time_delta_plus = TimeDelta(float(self.cfg_par['rfi']['chunks']['time_step'])*60., format='sec')
                    start = self.cfg_par['rfi']['startdate']+time_del
                    end = start+time_delta_plus
                    
                    self.logger.info((" ------ Plotting chunk #{0:d}:").format(i))
                    self.logger.info(("\t \t between {0:%d}{0:%b}{0:%y}: {0:%H}:{0:%M} - {1:%H}:{1:%M}").format(start.datetime,end.datetime))

                    if self.cfg_par['rfi']['rfi_enable']==False:

                        results = rfi.load_from_ms(self.cfg_par,timez)
                        self.logger.info("---- MSfile Loaded -----\n")    
                        if results != 1:
                            rfi.baselines_from_ms(self.cfg_par)
                            self.logger.info("---- Dataset sorted by baseline lenght ----\n")
                        else:
                            self.logger.warning("---- This chunk is empty ----\n")
                            continue       

                    rfiPL.plot_rfi_imshow(self.cfg_par,i)
                    self.logger.info("---- RFI in 2D plotted ----\n")
                    self.cfg_par['plots']['plot_noise'] = 'rfi'
                    self.cfg_par['plots']['long_short'] = False
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['long_short'] = True
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_noise'] = 'noise_factor'
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_noise'] = 'noise'
                    rfiPL.plot_noise_frequency(self.cfg_par,i)         
                    self.logger.info("---- RFI in 1D plotted ----\n")
                
                rfiPL.plot_altaz(self.cfg_par,68)
                self.logger.info("---- RFI in ALT/AZ plotted ----\n")
        
                if (self.cfg_par['plots']['movies']['altaz_gif']==True or self.cfg_par['plots']['movies']['2d_gif']==True or 
                    self.cfg_par['plots']['movies']['1d_gif']==True):
                    self.logger.info("---- Making movies ----\n")
     
                if self.cfg_par['plots']['movies']['altaz_gif']==True:

                    out_animation = self.cfg_par['general']['moviedir']+'AltAz_movie.gif'
                    filenames = rfi_files.find_altaz_plots(self.cfg_par)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                    self.logger.info("---- AltAz movie done ----\n")
                
                if self.cfg_par['plots']['movies']['2d_gif']==True:
                    out_animation = self.cfg_par['general']['moviedir']+'Time_2Dplot_movie.gif'
                    filenames = rfi_files.find_2d_plots(self.cfg_par)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                
                    self.logger.info("---- 2D movie done ----\n")
                
                if self.cfg_par['plots']['movies']['1d_gif']==True:
                    out_animation = self.cfg_par['general']['moviedir']+'TimeChunks_1D_flags.gif'
                    root_name = 'flags'
                    filenames = rfi_files.find_1d_plots(self.cfg_par,root_name)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                    
                    out_animation = self.cfg_par['general']['moviedir']+'TimeChunks_1D_noise.gif'
                    root_name = 'noise'
                    filenames = rfi_files.find_1d_plots(self.cfg_par,root_name)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)

                    out_animation = self.cfg_par['general']['moviedir']+'TimeChunks_1D_noisefactor.gif'
                    root_name = 'noisefactor'
                    filenames = rfi_files.find_1d_plots(self.cfg_par,root_name)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)         
                    self.logger.info("---- 1D movies done ----\n")
                
                if (self.cfg_par['plots']['movies']['altaz_gif']==True or self.cfg_par['plots']['movies']['2d_gif']==True or 
                    self.cfg_par['plots']['movies']['1d_gif']==True):
                    self.logger.info("---- Movies done ----\n")

                rfi_files.write_html_timereport(self.cfg_par)                 

            else:
                
                if self.cfg_par['rfi']['rfi_enable']==False:

                    results = rfi.load_from_ms(self.cfg_par,0)
                    self.logger.info("---- MSfile Loaded -----\n")    
                    if results != 1:
                        rfi.baselines_from_ms(self.cfg_par)
                        self.logger.info("---- Dataset sorted by baseline lenght ----\n")
                    else:
                        self.logger.warning("---- This dataset is empyy ----\n")

                rfiPL.plot_altaz_short(self.cfg_par)
                self.logger.info("---- Alt/Az plotted ----\n")                            
                rfiPL.plot_rfi_imshow(self.cfg_par,-1)
                self.logger.info("---- RFI in 2D plotted ----\n")
                self.cfg_par['plots']['plot_noise'] = 'rfi'
                self.cfg_par['plots']['long_short'] = False
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['long_short'] = True
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_noise'] = 'noise_factor'
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_noise'] = 'noise'
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.logger.info("---- RFI in 1D plotted ----\n")

                rfi_files.write_html_fullreport(self.cfg_par)

        self.logger.info(" ---- cleaning up ---- \n")

        rfi_files.cleanup(self.cfg_par)

        self.logger.info(" ---- End of RFInder ---- \n\n")

        return 0
