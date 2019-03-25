# Import modules
import os
import sys
import string
import numpy as np
import yaml
import json
import glob
import argparse
from  argparse import ArgumentParser
import textwrap as _textwrap
import logging
import logging.config

import warnings
import pkg_resources

from astropy.io import fits, ascii
from astropy import units as u
from astropy.time import Time, TimeDelta
from astropy.table import Table, Column, MaskedColumn



# get rfinder install directory
RFINDER_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RFINDER_DIR = RFINDER_PATH+'/rfinder/'

sys.path.append(os.path.join(RFINDER_PATH, 'rfinder'))

import rfi
import rfinder_stats as rfi_stats
import rfinder_plots as rfi_plots
import rfinder_files as rfiFL

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
    __version__ = pkg_resources.require("rfinder")[0].version
except pkg_resources.DistributionNotFound:
    __version__ = "dev"

####################################################################################################

class MultilineFormatter(argparse.HelpFormatter):
    def _fill_text(self, text, width, indent):
        text = self._whitespace_matcher.sub(' ', text).strip()
        paragraphs = text.split('|n ')
        multiline_text = ''
        for paragraph in paragraphs:
            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
            multiline_text = multiline_text + formatted_paragraph
        return multiline_text


def is_valid_file(parser, arg):
    if not os.path.exists(arg):
        parser.error("The file '%s' does not exist!" % arg)

    return arg



class rfinder:
    '''

    Class to investigate the RFI behaviour during observations

    '''

    def __init__(self):
        '''

        Set self.logger for spectrum extraction
        Find config file
        If not specified by user load rfinder_default.yml

        '''

        # set self.logger
        #with open(RFINDER_DIR+'/templates/logcfg.yml', 'r') as f:
        #    config = yaml.safe_load(f.read())
        #    logging.config.dictConfig(config)

        self.logger = logging.getLogger('log-rfinder.log')
        #self.logger.setLevel(logging.INFO)

        #fh = logging.FileHandler('log-rfinder.log')
        #fh.setLevel(logging.INFO)

        #ch = logging.StreamHandler()
        #ch.setLevel(logging.WARNING)

        #formatter = logging.Formatter('%(levelname)s - %(filename)s - %(message)s')
        #fh.setFormatter(formatter)
        #ch.setFormatter(formatter)

        #self.logger.addHandler(ch)
        #self.logger.addHandler(fh)


    def readArgs(self,argv):


        self.parser = ArgumentParser(description='RFInder: package to visualize the flagged RFI in a dataset '
                                '|n version {:s} |n install path {:s} |n '
                                'Filippo Maccagni <filippo.maccagni@gmial.com>'.format(__version__,
                                                                                   os.path.dirname(__file__)),
                                formatter_class=MultilineFormatter,
                                add_help=False)

        add = self.parser.add_argument

        add("-h", "--help",  action="store_true",
                help="Print help message and exit")

        add("-v","--version", action='version',
                version='{:s} version {:s}'.format(self.parser.prog, __version__))

        add('-c', '--config',
            type=lambda a: is_valid_file(self.parser, a),
            default=False,
            help='RFInder configuration file (YAML format)')

        add('-idir', '--input_dir',
            type= str,
            default = False,
            help= 'select working directory (MS file assumed to be here)')

        add('-odir', '--output_dir',
            type=str,
            default=False,
            help='select output directory')

        add('-l','--label',
            type=str,
            default=False,
            help='ouput folder is called: rfi_polarization_label')

        add('-i', '--input',
            type=str,
            default=False,
            help='''input ['MS'] file''')

        add('-fl', '--field',
            type=int,
            default=False,
            help='select field of MS file to analyze')

        add('-tel', '--telescope',
            type=str,
            default=False,
            help='select telescope: meerkat, apertif, wsrt')

        add('-mode', '--rfimode',
            type=str,
            default=False,
            help='select mode where to investigate RFI: use_flags or rms_clip')

        add('-pol', '--polarization',
            type=str,
            default=False,
            help='select stokes parameter: xx, yy, xy, yx, q (also in CAPS)')

        add('-fint', '--frequency_interval',
            nargs='*',
            default=False,
            help='select frequency interval where to measure noise in GHz')

        add('-spwAv', '--spw_av',
            type=int,
            default=False,
            help='select number of channels to average')

        add('-tStep', '--time_step',
            type=int,
            default=False,
            help='select time step in minutes in which divide the analysis of the MSfile')

        add('-sig', '--sigma_clip',
            type=int,
            default=False,
            help='select sigma clip for rms_clip mode to find RFI')

        add('-baseCut', '--baseline_cut',
            type=int,
            default=False,
            help='select cut in baseline lenght [m] for differential RFI analysis')

        add('-noCh', '--no_chunks',
            action='store_true',
            help='desable chunking in time')

        add('-yesCh','--yes_chunks',
            action='store_true',
            help='enable chunking in time')

        add('-noMov', '--no_movies',
            action='store_true',
            help='disable movies (use if dataset is read as a whole)')

        add('-noSpw', '--no_spw_av',
            action='store_true',
            help='disable averaging in channels')

        add('-yesSpw','--yes_spw_av',
            action='store_true',
            help='enable averaging in channels')

        add('-noClp','--no_cleanup',
            action='store_true',
            help='disable cleanup of intermediate products')

        add('-yesClp','--yes_cleanup',
            action='store_true',
            help='enable cleanup of intermediate products')



        args = self.parser.parse_args(argv)
        
        return args


    def setArgs(self,args):

        if args.input_dir:
            self.cfg_par['general']['workdir'] = args.input_dir
        if args.output_dir:
            self.cfg_par['general']['outdir'] = args.output_dir
        if args.input:
            self.cfg_par['general']['msname'] = args.input
        if args.field != None:
            self.cfg_par['general']['field'] = args.field
        if args.telescope:
            self.cfg_par['general']['telescope'] = args.telescope
        if args.polarization:
            self.cfg_par['rfi']['polarization'] = args.polarization
        if args.baseline_cut:
            self.cfg_par['rfi']['baseline_cut'] = args.baseline_cut
        if args.time_step:
            self.cfg_par['rfi']['chunks']['time_enable'] = True
            self.cfg_par['rfi']['chunks']['time_step'] = args.time_step
        if args.spw_av:
            self.cfg_par['rfi']['chunks']['spw_enable'] = True
            self.cfg_par['rfi']['chunks']['spw_width'] = args.spw_av
 
        if args.no_chunks==True:
                self.cfg_par['rfi']['chunks']['time_enable'] = False
        if args.yes_chunks == True:
                self.cfg_par['rfi']['chunks']['time_enable'] = True

        if args.no_spw_av==True:
            self.cfg_par['rfi']['chunks']['spw_enable'] = False
        
        if args.yes_spw_av==True:
            self.cfg_par['rfi']['chunks']['spw_enable'] = True

        if args.no_movies == True : 
            self.cfg_par['plots']['movies']['movies_in_report'] = False

        if args.no_cleanup==True:
            self.cfg_par['general']['cleanup_enable'] = False
        if args.yes_cleanup==True:
            self.cfg_par['general']['cleanup_enable'] = True
        
        if args.label !=False :
            self.cfg_par['general']['outlabel'] = '_'+args.label
        
        else:
            self.cfg_par['general']['outlabel'] = ''               
        
        if (args.rfimode == 'rms_clip' or args.rfimode == 'rms') :
                self.cfg_par['rfi']['RFInder_mode'] = args.rfimode
                if args.sigma_clip:
                    self.cfg_par['rfi']['rms_clip'] = args.sigma_clip
                if args.frequency_interval:
                    self.cfg_par['rfi']['noise_measure_edges'] = args.frequency_interval

        return self

    def set_cfg_par(self):

        key = 'general'

        workdir  = self.cfg_par[key].get('workdir', None)
        msfile = workdir + self.cfg_par[key].get('msname', None)
        self.cfg_par[key]['msfullpath'] = msfile      

        outdir  = self.cfg_par[key].get('outdir', None)
        rfidir  = outdir+'rfi_'+self.cfg_par['rfi']['polarization']+self.cfg_par[key]['outlabel']+'/'
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

        self.logger.warning("------ STARTING RFI analysis ------\n")

        if self.cfg_par[task]['rfi_enable']==True:
            
            if self.cfg_par[task]['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.warning("------ Working on time chunks ------\n")

                for i in xrange(0,len(times)-1):
                    timez = [times[i],times[i+1]] 
                    
                    #time chunk properties
                    time_delta = float(self.cfg_par['rfi']['chunks']['time_step'])*i
                    time_del = TimeDelta(time_delta*60., format='sec')
                    time_delta_plus = TimeDelta(float(self.cfg_par['rfi']['chunks']['time_step'])*60., format='sec')
                    start = self.cfg_par['rfi']['startdate']+time_del
                    end = start+time_delta_plus
                    self.logger.warning((" ------ Working on chunk #{0:d}:").format(i))
                    self.logger.warning(("\tbetween {0:%d}{0:%b}{0:%y}: {0:%H}:{0:%M} - {1:%H}:{1:%M}\n").format(start.datetime,end.datetime))

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
                        self.logger.warning(" ------ This chunk is empty ------\n")
                        continue

                self.logger.info("------ End of RFI analysis on time chunks ------\n")

            else:

                rfi.load_from_ms(self.cfg_par,0,0)
                #determine alt/az

                self.logger.warning("------ MSfile Loaded -----\n")
                rfi.baselines_from_ms(self.cfg_par)
                self.logger.warning("------ Dataset sorted by baseline lenght ------\n")
                datas = rfi.priors_flag(self.cfg_par)
                self.logger.warning("------ Bad antennas and autocorrelations flagged ------\n")
                rfi.find_rfi(datas,self.cfg_par,-1)
                self.logger.warning("------  RFI found  ------\n")
                rfiFL.rfi_frequency(self.cfg_par,-1)
                self.logger.warning("------ RFI saved to table ------\n")
                self.logger.warning("------ End of RFI analysis ------\n")
      
        task = 'plots'
        if self.cfg_par[task]['plot_enable']==True:
            
            if self.cfg_par['rfi']['chunks']['time_enable']==True:

                times, start, end = rfiST.time_chunk(self.cfg_par)
                self.logger.warning(" ------ Plotting on time chunks ------\n")

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

                        results = rfi.load_from_ms(self.cfg_par,timez,i)
                        self.logger.info("------ MSfile Loaded ------\n")    
                        if results != 1:
                            rfi.baselines_from_ms(self.cfg_par)
                            self.logger.info("------ Dataset sorted by baseline lenght ------\n")
                        else:
                            self.logger.warning("------ This chunk is empty ------\n")
                            continue       

                    rfiPL.plot_rfi_imshow(self.cfg_par,i)
                    self.logger.info("------ RFI in 2D plotted ------\n")
                    self.cfg_par['plots']['plot_noise'] = 'rfi'
                    self.cfg_par['plots']['long_short'] = False
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['long_short'] = True
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_noise'] = 'noise_factor'
                    rfiPL.plot_noise_frequency(self.cfg_par,i)
                    self.cfg_par['plots']['plot_noise'] = 'noise'
                    rfiPL.plot_noise_frequency(self.cfg_par,i)         
                    self.logger.warning("------ RFI in 1D plotted ------\n")
                
                rfiPL.plot_altaz(self.cfg_par,68)
                self.logger.warning("------ RFI in ALT/AZ plotted ------\n")
        
                if (self.cfg_par['plots']['movies']['altaz_gif']==True or self.cfg_par['plots']['movies']['2d_gif']==True or 
                    self.cfg_par['plots']['movies']['1d_gif']==True):
                    self.logger.warning("------ Making movies ------\n")
     
                if self.cfg_par['plots']['movies']['altaz_gif']==True:

                    out_animation = self.cfg_par['general']['moviedir']+'AltAz_movie.gif'
                    filenames = rfiFL.find_altaz_plots(self.cfg_par)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                    self.logger.info("------ AltAz movie done ------\n")
                
                if self.cfg_par['plots']['movies']['2d_gif']==True:
                    out_animation = self.cfg_par['general']['moviedir']+'Time_2Dplot_movie.gif'
                    filenames = rfiFL.find_2d_plots(self.cfg_par)
                    rfiPL.gif_me_up(self.cfg_par,filenames,out_animation)
                
                    self.logger.info("------ 2D movie done ------\n")
                
                if self.cfg_par['plots']['movies']['1d_gif']==True:
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
                
                if (self.cfg_par['plots']['movies']['altaz_gif']==True or self.cfg_par['plots']['movies']['2d_gif']==True or 
                    self.cfg_par['plots']['movies']['1d_gif']==True):
                    self.logger.warning("------ Movies done ------\n")

                rfiFL.write_html_timereport(self.cfg_par)                 

            else:
                
                if self.cfg_par['rfi']['rfi_enable']==False:

                    results = rfi.load_from_ms(self.cfg_par,0,0)
                    self.logger.info("------ MSfile Loaded -----\n")    
                    if results != 1:
                        rfi.baselines_from_ms(self.cfg_par)
                        self.logger.info("------ Dataset sorted by baseline lenght ------\n")
                    else:
                        self.logger.info("------ This dataset is empyy ------\n")

                rfiPL.plot_altaz_short(self.cfg_par)
                self.logger.warning("------ Alt/Az plotted ------\n")                            
                rfiPL.plot_rfi_imshow(self.cfg_par,-1)
                self.logger.warning("------ RFI in 2D plotted ------\n")
                self.cfg_par['plots']['plot_noise'] = 'rfi'
                self.cfg_par['plots']['long_short'] = False
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['long_short'] = True
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_noise'] = 'noise_factor'
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.cfg_par['plots']['plot_noise'] = 'noise'
                rfiPL.plot_noise_frequency(self.cfg_par,-1)
                self.logger.warning("------ RFI in 1D plotted ------\n")

                rfiFL.write_html_fullreport(self.cfg_par)

        self.logger.info("------ cleaning up ------\n")
        
        if self.cfg_par['general']['cleanup_enable'] == True:
            rfiFL.cleanup(self.cfg_par)

        self.logger.warning("------ End of RFInder ------\n\n")

        return 0



    def main (self,argv):
        
        for i, arg in enumerate(argv):
            if (arg[0] == '-') and arg[1].isdigit(): argv[i] = ' ' + arg

        args = self.readArgs(argv)

        if args.help:  #rfinder -h 
            self.parser.print_help()

            print("""\nRun a command. This can be:\n \nrfinder \nrfinder -c path_to_config_file.yml
rfinder -i <ngc1399.ms> -fl <num> -tel <meerkat/apertif/wsrt>\n""")

            sys.exit(0)

        elif args.config:    #rfinder -c config_file.yml         
            self.logger.warning('------ Reading your parameter file ------\n')
            # read database here
            files =  args.config
            cfg = open(files)
            self.cfg_par = yaml.load(cfg)

        else: #rfinder  or rfinder -options
            workdir = os.getcwd()
            workdir = workdir+'/'
            exists = os.path.isfile(workdir+'/'+DEFAULT_CONFIG)
            if exists:
                self.logger.warning('------ Reading default parameter file in your directory ------\n')
                file_default = os.path.join(workdir, DEFAULT_CONFIG)
                cfg = open(file_default)
                self.cfg_par = yaml.load(cfg)
            else:
                # Keep presets
                self.logger.warning('------ Reading default installation parameter file ------\n')
                file_default = os.path.join(RFINDER_DIR, DEFAULT_CONFIG)
                cfg = open(file_default)
                self.cfg_par = yaml.load(cfg)            
                self.cfg_par['general']['workdir'] = workdir
                self.cfg_par['general']['outdir'] = workdir
                with open(workdir+DEFAULT_CONFIG, 'w') as outfile:
                    yaml.dump(self.cfg_par, outfile, default_flow_style=False)

                if args.field ==False and args.input==False :

                    self.logger.warning('''MSNAME & telescope missing
              \t\tplease edit rfinder_default.yml in your current directory
              \t\tor run: rfinder -i msname -fl <field_number> -tel <meerkat,apertif,wsrt>
              \t\t(assuming the observation is located in your current directory)
                    \n''')
                    self.logger.critical('''------ RFInder out ------\n''')
                    
                    sys.exit(0)

                else:
                    self.logger.warning('''------ you provided MSname and telescope in your first run, 
            \tassuming MS is your current directory ------\n''')
 
            if all(x==False for x in vars(args).values()) == False and (args.help==False and args.config == False):
                self.logger.warning('------ Updating arguments given from terminal ------\n')

                self.setArgs(args)
                with open(workdir+DEFAULT_CONFIG, 'w') as outfile:
                    yaml.dump(self.cfg_par, outfile, default_flow_style=False)
 
        self.cfg_par['general']['template_folder'] = os.path.join(RFINDER_PATH,'rfinder/templates')
        self.set_cfg_par()

        return self
