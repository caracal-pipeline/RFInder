# Import modules
import os
import sys
import string
import numpy as np
import yaml
import json
import glob
import click
from  argparse import ArgumentParser
import textwrap as _textwrap
import logging
import logging.config

import warnings
from importlib.metadata import version, PackageNotFoundError

from astropy.io import fits, ascii
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.table import Table, Column, MaskedColumn

from scabha.schema_utils import clickify_parameters
from omegaconf import OmegaConf

# get rfinder install directory
RFINDER_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RFINDER_DIR = RFINDER_PATH+'/rfinder/'

sys.path.append(os.path.join(RFINDER_PATH, 'rfinder'))

from rfinder import rfi
from rfinder import rfinder_stats as rfi_stats
from rfinder import rfinder_plots as rfi_plots
from rfinder import rfinder_files as rfiFL


schemas = OmegaConf.load(os.path.join(os.path.dirname(__file__), "rfinder.yaml"))

rfi = rfi.rfi()
rfiST = rfi_stats.rfi_stats()
rfiPL = rfi_plots.rfi_plots()

DEFAULT_CONFIG = 'rfinder_default.yml'
# Set up logging infrastructure
LOG_FILE = 'log-rfinder.log'
if os.path.exists(LOG_FILE) == True:
    os.remove(LOG_FILE)

# This is is the default log file. It logs stimela images, containers and processes

if not sys.warnoptions:
    warnings.simplefilter("ignore")

try:
        __version__ = version("rfinder")
except PackageNotFoundError:
        __version__ = "dev"


####################################################################################################


class Rfinder:
    '''

    Class to investigate the RFI behaviour during observations

    '''

    def __init__(self):
        '''

        Set self.logger for spectrum extraction
        Find config file
        If not specified by user load rfinder_default.yml

        '''

        self.logger = logging.getLogger('log-rfinder.log')
        self.logger.setLevel(logging.INFO)


    def setArgs(self, kwargs):
        if kwargs.get('input_dir'):
            self.cfg_par['general']['workdir'] = kwargs['input_dir']
        if kwargs.get('output_dir'):
            self.cfg_par['general']['outdir'] = kwargs['output_dir']
        if kwargs.get('msname'):
            self.cfg_par['general']['msname'] = kwargs['msname']
        if kwargs.get('field'):
            self.cfg_par['general']['field'] = kwargs['field']
        if kwargs.get('ncpu'):
            self.cfg_par['general']['ncpu'] = kwargs['ncpu']
        if kwargs.get('telescope'):
            self.cfg_par['general']['telescope']['name'] = kwargs['telescope']
        if kwargs.get('polarization'):
            self.cfg_par['rfi']['polarization'] = kwargs['polarization']
        if kwargs.get('baseline_cut'):
            self.cfg_par['rfi']['baseline_cut'] = kwargs['baseline_cut']
        if kwargs.get('chunks_time_step'):
            self.cfg_par['rfi']['chunks']['time_enable'] = True
            self.cfg_par['rfi']['chunks']['time_step'] = kwargs['chunks_time_step']
        if kwargs.get('chunks_spw_width'):
            self.cfg_par['rfi']['chunks']['spw_enable'] = True
            self.cfg_par['rfi']['chunks']['spw_width'] = kwargs['chunks_spw_width']

        if kwargs.get('chunks_time_enable'):
            self.cfg_par['rfi']['chunks']['time_enable'] = True
        elif kwargs.get('chunks_time_enable') is False:
            self.cfg_par['rfi']['chunks']['time_enable'] = False

        if kwargs.get('chunks_spw_enable') is True:
            self.cfg_par['rfi']['chunks']['spw_enable'] = True
        elif kwargs.get('chunks_spw_enable') is False:
            self.cfg_par['rfi']['chunks']['spw_enable'] = False

        if kwargs.get('plot_details_movies_in_report') is True:
            self.cfg_par['plots']['plot_details']['movies']['movies_in_report'] = True
        elif kwargs.get('plot_details_movies_in_report') is False:
            self.cfg_par['plots']['plot_details']['movies']['movies_in_report'] = False

        if kwargs.get('cleanup_enable') is True:
            self.cfg_par['general']['cleanup_enable'] = True
        elif kwargs.get('cleanup_enable') is False:
            self.cfg_par['general']['cleanup_enable'] = False

        if kwargs.get('label'):
            self.cfg_par['general']['outlabel'] = '_' + kwargs['label']
        else:
            self.cfg_par['general']['outlabel'] = '_' + self.cfg_par['general']['outlabel']
        if kwargs.get('rfimode') in ['rms_clip', 'use_flags']:
            self.cfg_par['rfi']['RFInder_mode'] = kwargs['rfimode']
            self.cfg_par['rfi']['rfi_enable'] = True
            if kwargs.get('rms_clip'):
                self.cfg_par['rfi']['rms_clip'] = kwargs['sigma_clip']
            if kwargs.get('frequency_interval'):
                self.cfg_par['rfi']['noise_measure_edges'] = kwargs['frequency_interval']
            if kwargs.get('plot_details_enable'):
                self.cfg_par['plots']['plot_details']['enable'] = True
        if kwargs.get('rfi_enable') is True:
            self.cfg_par['rfi']['rfi_enable'] = True
            self.cfg_par['plots']['plot_details']['enable'] = True
        elif kwargs.get('rfi_enable') is False:
            self.cfg_par['rfi']['rfi_enable'] = False
            self.cfg_par['plots']['plot_details']['enable'] = False
        if kwargs.get('plot_summary_enable') is True:
            self.cfg_par['plots']['plot_summary']['enable'] = True
            if kwargs.get('plot_summary_axis'):
                self.cfg_par['plots']['plot_summary']['axis'] = kwargs['plot_summary_axis']
            if kwargs.get('plot_summary_freq_bin'):
                self.cfg_par['plots']['plot_summary']['freq_bin'] = kwargs['plot_summary_freq_bin']
            if kwargs.get('plot_summary_report'):
                self.cfg_par['plots']['plot_summary']['report'] = kwargs['plot_summary_report']
        elif kwargs.get('plot_summary_enable') in [False, None]:
            if kwargs.get('plot_summary_enable') is False:
                self.cfg_par['plots']['plot_summary']['enable'] = False
            if kwargs.get('plot_summary_axis'):
                self.cfg_par['plots']['plot_summary']['axis'] = kwargs['plot_summary_axis']
            if kwargs.get('plot_summary_freq_bin'):
                self.cfg_par['plots']['plot_summary']['freq_bin'] = kwargs['plot_summary_freq_bin']
            if kwargs.get('plot_summary_report'):
                self.cfg_par['plots']['plot_summary']['report'] = kwargs['plot_summary_report']

        return self

    def set_cfg_par(self):

        key = 'general'

        workdir  = self.cfg_par[key].get('workdir', None)
        msfile = workdir + self.cfg_par[key].get('msname', None)
        self.cfg_par[key]['msfullpath'] = msfile      

        outdir  = self.cfg_par[key].get('outdir', None)
        rfidir  = outdir+'rfi_'+self.cfg_par['rfi']['polarization']+'_'+self.cfg_par[key]['outlabel']+'/'
        self.cfg_par[key]['rfidir'] = rfidir

        tabledir = rfidir+'tables/'
        self.cfg_par[key]['tabledir'] = tabledir
    
        rfiplotdir = rfidir+'plots/'
        self.cfg_par[key]['plotdir'] = rfiplotdir 

        moviedir = rfidir+'plots/movies/'
        self.cfg_par[key]['moviedir'] = moviedir             

        rfitimedir = rfidir+'time_chunks/'
        self.cfg_par[key]['rfitimedir'] = rfitimedir

        timetabledir = tabledir+'time_chunks/'
        self.cfg_par[key]['timetabledir'] = timetabledir

        timeplotdir_tmp = rfiplotdir+'time_chunks/'
        self.cfg_par[key]['timechunksdir'] = timeplotdir_tmp          

        timeplotdir1d = timeplotdir_tmp+'1D/'
        self.cfg_par[key]['timeplotdir1D'] = timeplotdir1d

        timeplotdir2d = timeplotdir_tmp+'2D/'
        self.cfg_par[key]['timeplotdir2D'] = timeplotdir2d

        altazplotdir = rfidir+'plots/altaz/'
        self.cfg_par[key]['altazplotdir'] = altazplotdir        


        return 0




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
            4: 1d plot of overall RFI flagged by scan, antenna and correlation.
        If cfg_par['rfi']['chunks']['time_chunks'] is enabled
            1: executes 'rfi' and 'plots' procedure dividing the observation in time-steps given by cfg_par['rfi']['chunks']['time_step']
            2: collects the info about the % of RFI for each time step in Alt/Az plots 
        If cfg_par['beam_shape'] is enabled
            1: create FLAG column in MS file (rfi_flag)
            2: determine psf using wsclean (make_psf)
        '''

        # cont_sources

        self.cfg_par = cfg_par
        task = 'rfi'
        rfiFL.set_dirs(self.cfg_par)

        self.logger.info("------ STARTING RFI analysis ------\n")

        if self.cfg_par[task]['rfi_enable']==True:
            
            if self.cfg_par[task]['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.info("------ Working on time chunks ------\n")

                for i in range(len(times)-1):
                    timez = [times[i],times[i+1]] 
                    
                    #time chunk properties
                    time_delta = float(self.cfg_par['rfi']['chunks']['time_step'])*i
                    time_del = TimeDelta(time_delta*60., format='sec')
                    time_delta_plus = TimeDelta(float(self.cfg_par['rfi']['chunks']['time_step'])*60., format='sec')
                    start = self.cfg_par['rfi']['startdate']+time_del
                    end = start+time_delta_plus
                    self.logger.info((" ------ Working on chunk #{0:d}:").format(i))
                    self.logger.info(("\tbetween {0:%d}{0:%b}{0:%y}: {0:%H}:{0:%M} - {1:%H}:{1:%M}\n").format(start.datetime,end.datetime))

                    result = rfi.load_from_ms(self.cfg_par,timez,i)
                    self.logger.info("------ MSfile Loaded ------\n")

                    #sort visibilities by baseline lenght
                    if result != 1:
                        rfi.baselines_from_ms(self.cfg_par)
                        self.logger.info("------ Dataset sorted by baseline lenght ------\n")

                        #flag bad antennas (from configuration file)
                        datas = rfi.priors_flag(self.cfg_par)
                        self.logger.info("------ Bad antennas and autocorrelations flagged ------\n")

                        #find rfi above threshold
                        rfi.find_rfi(datas,self.cfg_par,i)
                        self.logger.info(" ------  RFI found  ------\n")

                        rfiFL.rfi_frequency(self.cfg_par,i)
                        self.logger.info("------ RFI saved to table ------\n")
                    else:
                        self.logger.info(" ------ This chunk is empty ------\n")
                        continue

                self.logger.info("------ End of RFI analysis on time chunks ------\n")

            else:

                rfi.load_from_ms(self.cfg_par,0,0)
                #determine alt/az

                self.logger.info("------ MSfile Loaded -----\n")
                rfi.baselines_from_ms(self.cfg_par)
                self.logger.info("------ Dataset sorted by baseline lenght ------\n")
                datas = rfi.priors_flag(self.cfg_par)
                self.logger.info("------ Bad antennas and autocorrelations flagged ------\n")
                rfi.find_rfi(datas,self.cfg_par,-1)
                self.logger.info("------  RFI found  ------\n")
                rfiFL.rfi_frequency(self.cfg_par,-1)
                self.logger.info("------ RFI saved to table ------\n")
                self.logger.info("------ End of RFI analysis ------\n")

        else:
            rfi.load_from_ms(self.cfg_par,0,0)
            self.logger.info("------ MSfile Loaded -----\n")
      
        task = 'plots'

        if self.cfg_par[task]['plot_details']['enable']==True:

            if self.cfg_par['rfi']['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.info(" ------ Plotting on time chunks ------\n")

                for i in range(len(times)-1):

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

                        results = rfi.load_from_ms(self.cfg_par,timez,i)
                        self.logger.info("------ MSfile Loaded ------\n")    
                        if results != 1:
                            rfi.baselines_from_ms(self.cfg_par)
                            self.logger.info("------ Dataset sorted by baseline lenght ------\n")
                        else:
                            self.logger.info("------ This chunk is empty ------\n")
                            continue       

                    rfiPL.plot_rfi_imshow(self.cfg_par,i)
                    self.logger.info("------ RFI in 2D plotted ------\n")
                    self.cfg_par['plots']['plot_details']['plot_noise'] = 'rfi'
                    self.cfg_par['plots']['plot_details']['long_short'] = False
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_details']['long_short'] = True
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_details']['plot_noise'] = 'noise_factor'
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_details']['plot_noise'] = 'noise'
                    rfiPL.plot_noise_frequency(self.cfg_par,i)         
                    self.logger.info("------ RFI in 1D plotted ------\n")
                
                rfiPL.plot_altaz(self.cfg_par,68)
                self.logger.info("------ RFI in ALT/AZ plotted ------\n")
        
                if (self.cfg_par['plots']['plot_details']['movies']['altaz_gif']==True or
                    self.cfg_par['plots']['plot_details']['movies']['2d_gif']==True or
                    self.cfg_par['plots']['plot_details']['movies']['1d_gif']==True):
                    self.logger.info("------ Making movies ------\n")
     
                if self.cfg_par['plots']['plot_details']['movies']['altaz_gif']==True:

                    out_animation = self.cfg_par['general']['moviedir']+'AltAz_movie.gif'
                    filenames = rfiFL.find_altaz_plots(self.cfg_par)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                    self.logger.info("------ AltAz movie done ------\n")
                
                if self.cfg_par['plots']['plot_details']['movies']['2d_gif']==True:
                    out_animation = self.cfg_par['general']['moviedir']+'Time_2Dplot_movie.gif'
                    filenames = rfiFL.find_2d_plots(self.cfg_par)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                
                    self.logger.info("------ 2D movie done ------\n")
                
                if self.cfg_par['plots']['plot_details']['movies']['1d_gif']==True:
                    out_animation = self.cfg_par['general']['moviedir']+'TimeChunks_1D_flags.gif'
                    root_name = 'flags'
                    filenames = rfiFL.find_1d_plots(self.cfg_par,root_name)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                    
                    out_animation = self.cfg_par['general']['moviedir']+'TimeChunks_1D_noise.gif'
                    root_name = 'noise'
                    filenames = rfiFL.find_1d_plots(self.cfg_par,root_name)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)

                    out_animation = self.cfg_par['general']['moviedir']+'TimeChunks_1D_noisefactor.gif'
                    root_name = 'noisefactor'
                    filenames = rfiFL.find_1d_plots(self.cfg_par,root_name)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)         
                    self.logger.info("------ 1D movies done ------\n")
                
                if (self.cfg_par['plots']['plot_details']['movies']['altaz_gif']==True or
                    self.cfg_par['plots']['plot_details']['movies']['2d_gif']==True or
                    self.cfg_par['plots']['plot_details']['movies']['1d_gif']==True):
                    self.logger.info("------ Movies done ------\n")

                rfiFL.write_html_timereport(self.cfg_par)                 

            else:
                
                if self.cfg_par['rfi']['rfi_enable']==False:

                    results = rfi.load_from_ms(self.cfg_par,0,0)
                    self.logger.info("------ MSfile Loaded -----\n")    
                    if results != 1:
                        rfi.baselines_from_ms(self.cfg_par)
                        self.logger.info("------ Dataset sorted by baseline lenght ------\n")
                    else:
                        self.logger.info("------ This dataset is empty ------\n")

                rfiPL.plot_altaz_short(self.cfg_par)
                self.logger.info("------ Alt/Az plotted ------\n")                            
                rfiPL.plot_rfi_imshow(self.cfg_par,-1)
                self.logger.info("------ RFI in 2D plotted ------\n")
                self.cfg_par['plots']['plot_details']['plot_noise'] = 'rfi'
                self.cfg_par['plots']['plot_details']['long_short'] = False
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_details']['long_short'] = True
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_details']['plot_noise'] = 'noise_factor'
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_details']['plot_noise'] = 'noise'
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.logger.info("------ RFI in 1D plotted ------\n")
                rfiFL.write_html_fullreport(self.cfg_par)

        if self.cfg_par[task]['plot_summary']['enable']==True:
            summary_results = {}

            for axis in  self.cfg_par[task]['plot_summary']['axis']:
                flag_stats = rfiST.get_flags_summary_stats(self.cfg_par, axis)
                summary_results[axis] = dict(flag_stats)
                self.logger.info(f" ------ Plotting {axis} summary plots ------\n")
                rfiPL.plot_summary_stats(flag_stats, self.cfg_par, axis)
                self.logger.info("------ Summary plot done ------\n")

            if summary_results:
                self.logger.info(f'------ Total % Flagged: {round(sum(summary_results[axis].values())/len(summary_results[axis].values()),2)} ------')
                json_file = cfg_par['general']['rfidir'] + f'summary.json'
                with open(json_file, 'w') as f:
                    json.dump(summary_results, f)

            rfiFL.write_html_summaryreport(self.cfg_par)

        self.logger.info("------ cleaning up ------\n")
        
        if self.cfg_par['general']['cleanup_enable'] == True:
            rfiFL.cleanup(self.cfg_par)

        self.logger.info("------ End of RFInder ------\n\n")

        return 0


    def main (self, **kwargs):
        config = kwargs.get('_config')

        if config:    #rfinder -c config_file.yml
            self.logger.info('------ Reading your parameter file ------\n')
            # read database here
            cfg = open(config)
            self.cfg_par = yaml.load(cfg, Loader=yaml.Loader)

        else: #rfinder  or rfinder -options
            workdir = os.getcwd()
            workdir = workdir+'/'
            exists = os.path.isfile(workdir+'/'+DEFAULT_CONFIG)
            if exists:
                self.logger.info('------ Reading default parameter file in your directory ------\n')
                file_default = os.path.join(workdir, DEFAULT_CONFIG)
                cfg = open(file_default)
                self.cfg_par = yaml.load(cfg, Loader=yaml.Loader)
            else:
                # Keep presets
                self.logger.info('------ Reading default installation parameter file ------\n')
                file_default = os.path.join(RFINDER_DIR, DEFAULT_CONFIG)
                cfg = open(file_default)
                self.cfg_par = yaml.load(cfg, Loader=yaml.Loader)
                self.cfg_par['general']['workdir'] = workdir
                self.cfg_par['general']['outdir'] = workdir
                with open(workdir+DEFAULT_CONFIG, 'w') as outfile:
                    yaml.dump(self.cfg_par, outfile, default_flow_style=False)

                if kwargs['field'] == None and kwargs['msname'] == None :

                    self.logger.warning('''MSNAME & telescope missing
              \t\tplease edit rfinder_default.yml in your current directory
              \t\tor run: rfinder -i msname -fl <field_number> -tel <meerkat,apertif,wsrt>
              \t\t(assuming the observation is located in your current directory)
                    \n''')
                    self.logger.critical('''------ RFInder out ------\n''')
                    
                    sys.exit(0)

                else:
                    self.logger.info('''------ you provided MSname and telescope in your first run, 
                                        \tassuming MS is your current directory ------\n''')
            if any(kwargs.values()) and not (kwargs.get('help') or kwargs.get('config')):
                self.logger.info('------ Updating arguments given from terminal ------\n')

                self.setArgs(kwargs)
                with open(workdir+DEFAULT_CONFIG, 'w') as outfile:
                    yaml.dump(self.cfg_par, outfile, default_flow_style=False)
 
        self.cfg_par['general']['template_folder'] = os.path.join(RFINDER_PATH,'rfinder/templates')
        self.set_cfg_par()

        return self


@click.command("rfinder")
@clickify_parameters(schemas.cabs.get("rfinder"))
@click.option('-h', '--help', is_flag=True, help="Show this message and exit.")
@click.version_option(version=__version__)
def driver(help, **kw):

    if help:  #rfinder -h
        print('RFInder: package to visualize the flagged RFI in a dataset\n'
              'version {:s}\n'
              'install path {:s}\n'
              'Filippo Maccagni <filippo.maccagni@gmial.com>\n'.format(__version__, os.path.dirname(__file__)))
        click.echo(driver.get_help(click.Context(driver)))
        print("\nRun a command. This can be:\n \nrfinder \nrfinder -c path_to_config_file.yml" +
                "\nrfinder -i <ngc1399.ms> -fl <num> -tel <meerkat/apertif/wsrt>" +
                "\nrfinder -i <ngc1399.ms> -fl <num> -tel <meerkat/apertif/wsrt> -rfi -mode rms_clip")
        sys.exit(0)

    logger = logging.getLogger('log-rfinder.log')
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.debug('info')
    logger.info('info')

    log_dir = kw.get('output_dir') or '.'
    os.makedirs(log_dir, exist_ok=True)
    fh = logging.FileHandler(os.path.join(log_dir, LOG_FILE))
    fh.setLevel(logging.INFO)

    ch = logging.StreamHandler()
    ch.setLevel(logging.INFO)

    formatter = logging.Formatter('%(asctime)s; %(levelname)s - %(filename)s - %(message)s')
    formatter_ch = logging.Formatter('%(asctime)s; %(message)s')

    fh.setFormatter(formatter)
    ch.setFormatter(formatter_ch)

    logger.addHandler(ch)
    logger.addHandler(fh)


    RFInder = Rfinder()
    rfi_par = RFInder.main(**kw)

    run = rfi_par.go(rfi_par.cfg_par)

    if run == 0:
        logger.info('\t+------+\n\t  Done\n\t+------+')
