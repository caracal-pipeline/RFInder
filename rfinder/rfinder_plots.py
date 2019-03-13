import os,string,sys, glob
import numpy as np

import matplotlib
matplotlib.rcParams['text.usetex'] = True

from astropy import units as u

from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import ticker, rc
from matplotlib.ticker import NullFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

from IPython.display import HTML, display

import matplotlib.image as mgimg

from astropy.io import fits as fits
from astropy.table import Table
from astropy.time import Time, TimeDelta

import imageio as io


import logging

import rfi
rfi = rfi.rfi()
import rfinder_stats as rfi_stats
rfiST = rfi_stats.rfi_stats()


class rfi_plots:

    def __init__(self):

        #set self.logger
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

    # Project spherical coordinates on XY plane
    def aitoff(self,x,y):
        x=x/180*np.pi
        y=y/180*np.pi
        a=np.arccos(np.cos(y)*np.cos(x/2))
        x=2*np.cos(y)*np.sin(x/2)/np.sin(a)*a/np.pi*180
        y=np.sin(y)/np.sin(a)*a/np.pi*180
        
        return x,y

    def plot_rfi_imshow(self,cfg_par,time_step=-1):      
        '''
        
        Plots the .fits image output of rfi_im jpg format.
        
        '''        
        
        #check if image exists
        #plot image

        self.logger.info("\t ... Plotting RFI in 2D ... \n")


        if time_step != -1:
            outputdir = cfg_par['general']['rfitimedir']
            plotdir = cfg_par['general']['timeplotdir2D']
            time_delta = float(cfg_par['rfi']['chunks']['time_step'])*time_step
            time_tmp = int(time_delta)
            if time_tmp == 0:
                time_name = '00'+str(time_tmp)+'m'
            elif time_tmp <100:
                time_name = '0'+str(time_tmp)+'m'
            else:
                time_name= str(time_tmp)+'m'    
        else:
            outputdir = cfg_par['general']['rfidir']
            plotdir = cfg_par['general']['plotdir']        
            time_name = 'full' 

        rfi_freq_base = cfg_par['general']['rfidir']+'freq_base_'+time_name+'_im.fits'   

        if cfg_par['rfi']['RFInder_mode']== 'use_flags':
            rfi_freq_base =outputdir+'flags_base_'+time_name+'.fits'
            rfi_freq_base_plot = plotdir+'flags_base_'+time_name+'.png'
        if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
            rfi_freq_base =outputdir+'rfi_base_'+time_name+'.fits'
            rfi_freq_base_plot = plotdir+'rfi_base_'+time_name+'.png'
        
        #open file
        if os.path.exists(rfi_freq_base) == False:
            self.logger.error('### 2D RFI in fits format does not exist ###')    
        else:
            # open file
            t = fits.open(rfi_freq_base)
            data = t[0].data
            prihdr= t[0].header
            freqs = (np.linspace(1, data.shape[1], data.shape[1])- prihdr['CRPIX1'])*prihdr['CDELT1'] + prihdr['CRVAL1']
            
            # set x-y axes
            bandwidth= (cfg_par['rfi']['highfreq'] - cfg_par['rfi']['lowfreq'] )/1e6
            if bandwidth <=25:
                step_bin = 5
            elif bandwidth <=50 and bandwidth >25:
                step_bin = 10
            elif bandwidth <= 150. and bandwidth > 50.:
                step_bin=25.
            elif bandwidth <= 500. and bandwidth >  150.:
                step_bin=50.                
            else:
                step_bin=100.
            # Smaller multiple 
            freq_max = (freqs[-1] // step_bin) * step_bin  
            # Larger multiple 
            freq_min= (freqs[0] // step_bin) * step_bin + step_bin
            freqs_plot_bin=np.arange(freq_min,freq_max+step_bin,step_bin)
            freqs_plot_bin=np.round(freqs_plot_bin,0)
            freqs_plot_bin=np.array(freqs_plot_bin,dtype=int)

            freqs_plot_idx=np.zeros(freqs_plot_bin.shape)
            for i in xrange(0,len(freqs_plot_bin)):
                idx = (np.abs(freqs_plot_bin[i]-freqs)).argmin()
                freqs_plot_idx[i]=idx

            tele= cfg_par['general']['telescope']
            if tele == 'meerkat' or tele == 'MeerKAT' or tele == 'meerKAT' or tele == 'meer':    
                input_baselines = np.zeros(6)
                for i in xrange (0,len(input_baselines)):
                    input_baselines[i] = 100.*np.power(2,i)
            elif tele == 'apertif' or tele == 'Apertif' or tele == 'APERTIF' or tele == 'wsrt':
                input_baselines = np.zeros(7)
                for i in xrange (0,len(input_baselines)):
                    input_baselines[i] = 50.*np.power(2,i)


            input_baselines_idx=np.zeros(input_baselines.shape)

            baselines = np.array([cfg_par['rfi']['baseline_lenghts']])+0.
            
            if input_baselines[-1] > baselines[0,-1]:
                input_baselines[-1]=np.round(baselines[0,-1],0)

            for i in xrange(0, len(input_baselines)):
                idx = (np.abs(input_baselines[i] - baselines[0,:])).argmin()
                input_baselines_idx[i]=idx
            
            #set rc parameters 
            plt.rcParams['image.interpolation']='nearest'
            plt.rcParams['image.origin']='lower'
            plt.rcParams['image.aspect']='auto'

            params = {'font.family'         :' serif',
                      'font.style'          : 'normal',
                      'font.weight'         : 'book',
                      'font.size'           : 18.0,
                      'axes.linewidth'      : 1,
                      'lines.linewidth'     : 1,
                      'xtick.labelsize'     : 16,
                      'ytick.labelsize'     : 16, 
                      'xtick.direction'     :'in',
                      'ytick.direction'     :'in',
                      #'xtick.top'           : True,   # draw ticks on the top side
                      #'xtick.bottom'        : True,   # draw ticks on the bottom side    
                      #'ytick.left'          : True,   # draw ticks on the top side
                      #'ytick.right'         : True,   # draw ticks on the bottom side  
                      'xtick.major.size'    : 4,
                      'xtick.major.width'   : 1,
                      'xtick.minor.size'    : 2,
                      'xtick.minor.width'   : 1,
                      'ytick.major.size'    : 4,
                      'ytick.major.width'   : 1,
                      'ytick.minor.size'    : 2,
                      'ytick.minor.width'   : 1, 
                      'text.usetex'         : True,
                      'text.latex.unicode'  : True
                       }
            plt.rcParams.update(params)
            # Format axes
            nullfmt        = NullFormatter() 
            left, width    = 0.05, 0.9                                                                     #|These determine where the subplots go
            bottom, height = 0.1, 0.85
            #left_h = left+width+0.015
            #bottom_h = left+height+0.015

            box_centre    = [left, bottom, width, height]
            #box_x      = [left, bottom_h, width, 0.12]
            #box_y      = [left_h, bottom, 0.12, height]
            #box_cbar       = [left_h+0.13, bottom, 0.05, height]

            fig = plt.figure(1, figsize=(12,8))

            ax = plt.axes(box_centre)
            
            if bandwidth <= 500.:
                colormap = 'nipy_spectral_r'
            else:
                colormap = 'afmhot'


            im = ax.imshow(data,cmap=colormap,vmin=0,vmax=100)

            # ticks & labels
            ax.set_xticks(freqs_plot_idx)
            ax.set_yticks(input_baselines_idx)

            ax.set_xticklabels(freqs_plot_bin)
            ax.set_yticklabels(input_baselines)
                    
            ax.set_xlabel('Frequency [MHz]')
            ax.set_ylabel('Baseline lenght [m]')
            
            # colorbar

            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="3%", pad=0.05)
            cbar = fig.colorbar(im,ticks=[0,10,20, 30,40,50,60,70,80,90,100], cax=cax ) 

            
            if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
                cbar.set_label(r'$\% > 5 \times$ r.m.s.')
            if cfg_par['rfi']['RFInder_mode']== 'use_flags':
                cbar.set_label(r'$\%$ flagged visibilites')

            #times, start, end = rfi.time_chunk(cfg_par)
            #start.format='iso' 
            #start.subformat='date_hm'   
            #end.format='iso'
            #end.subformat='date_hm'


            if cfg_par['rfi']['chunks']['time_enable']== False:
                start = cfg_par['rfi']['startdate']
                end = cfg_par['rfi']['enddate']

            elif cfg_par['rfi']['chunks']['time_enable'] == True:
                time_del = TimeDelta(time_delta*60., format='sec')
                time_delta_plus = TimeDelta(float(cfg_par['rfi']['chunks']['time_step'])*60., format='sec')
                start = cfg_par['rfi']['startdate']+time_del
                end = start+time_delta_plus
           
            pol = str(cfg_par['rfi']['polarization'] )

            if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
                rfi_clip = str(cfg_par['rfi']['rms_clip'])+r'$\sigma$ clip'
                title_plot = 'Stokes:  {0:s} / {1:s} / {2:%d}{2:%b}{2:%y}: {2:%H}:{2:%M} - {3:%H}:{3:%M}'.format(pol,rfi_clip,start.datetime,end.datetime)
            if cfg_par['rfi']['RFInder_mode']== 'use_flags':
                title_plot = 'Stokes:  {0:s} / {1:s} / {2:%d}{2:%b}{2:%y}: {2:%H}:{2:%M} - {3:%H}:{3:%M}'.format(pol,'Flags',start.datetime,end.datetime)
            
            ax.set_title(title_plot)
         
            plt.savefig(rfi_freq_base_plot,format='png' ,overwrite=True)
            plt.close(fig)        
            self.logger.info("\t ... RFI per baseline lenght and frequency plotted ... \n\n")

    def plot_noise_frequency(self,cfg_par,time_step=-1):
        '''
        Plots the noise or % of rfi per frequency channel for all, long and short baselines.
        In default.cfga
            aperfi_noise = 'rfi' (or flag)
            aperfi_plot_long_short = False
        '''

        #open file
        #if os.path.exists(self.rfi_table) == False:
        #    self.self.logger.error('### Table of RFI and flags of visibilities does not exist ###')    
        #    self.self.logger.error('### Run aperfi.rfi_frequency() first ###')  
        #else:  

        self.logger.info("\t ... Plotting RFI in 1D ... \n")


        table_tmp = string.split(cfg_par['general']['msname'][0],'.MS')
        if len(table_tmp) == 1:
            table_tmp = string.split(cfg_par['general']['msname'][0],'.ms')

        if time_step != -1:
            tabledir = cfg_par['general']['timetabledir'] 
            plotdir = cfg_par['general']['timeplotdir1D']
            time_delta = float(cfg_par['rfi']['chunks']['time_step'])*time_step
            time_tmp = int(time_delta)
            #name file according to WSRT beams
            if time_tmp == 0:
                time_name = '00'+str(time_tmp)+'m'
            elif time_tmp <100:
                time_name = '0'+str(time_tmp)+'m'
            else:
                time_name= str(time_tmp)+'m'    
        else:
            tabledir = cfg_par['general']['tabledir']        
            plotdir = cfg_par['general']['plotdir']        
            time_name = 'full' 

        if cfg_par['rfi']['RFInder_mode']== 'use_flags':
            table_name = str(table_tmp[0])+'_flags_'+time_name+'.fits'
        if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
            table_name = str(table_tmp[0])+'_rfi_'+time_name+'.fits'

        rfi_table = tabledir+table_name
        
        #open file
        if os.path.exists(rfi_table) == False:
            self.logger.error('### Table of RFI results does not exist ###')    
        else:    
            
            t = fits.open(rfi_table)
            data_vec = t[1].data
            cols = t[1].columns
            
            freqs = np.array(data_vec['frequency'],dtype=float)
            flags = np.array(data_vec['percentage_flags'],dtype=float)
            noise_factor = np.array(data_vec['noise_factor'],dtype=float)
            noise_factor_long = np.array(data_vec['noise_factor_long'],dtype=float)
            flags_long = np.array(data_vec['percentage_flags_long'],dtype=float)
            noise_factor_short = np.array(data_vec['noise_factor_short'],dtype=float)
            flags_short = np.array(data_vec['percentage_flags_short'],dtype=float)

           
            if cfg_par['plots']['plot_noise'] == 'noise':
                rms = np.array(cfg_par['rfi']['theo_rms']*1e3,dtype=float)
                noise_all = noise_factor*rms
                noise_short = noise_factor_short*rms
                noise_long = noise_factor_long*rms
                out_plot = plotdir+'noise_'+time_name
            
            if cfg_par['plots']['plot_noise'] == 'noise_factor':
                self.logger.info("\t ... Plotting factor of noise increas per frequency channel ...")
                noise_all = noise_factor
                noise_short = noise_factor_short
                noise_long = noise_factor_long    
                out_plot = plotdir+'noisefactor_'+time_name
            
            if cfg_par['plots']['plot_noise'] == 'rfi':
                self.logger.info("\t ... Plotting percentage of flagged RFI per frequency channel ...")
                noise_all = flags
                noise_long = flags_long
                noise_short = flags_short
                out_plot = plotdir+'flags_'+time_name


                    # initialize plotting parameters
            params = {'font.family'         :' serif',
                      'font.style'          : 'normal',
                      'font.weight'         : 'book',
                      'font.size'           : 18.0,
                      'axes.linewidth'      : 1,
                      'lines.linewidth'     : 1,
                      'xtick.labelsize'     : 16,
                      'ytick.labelsize'     : 16, 
                      'xtick.direction'     :'in',
                      'ytick.direction'     :'in',
                      #'xtick.top'           : True,   # draw ticks on the top side
                      #'xtick.bottom'        : True,   # draw ticks on the bottom side    
                      #'ytick.left'          : True,   # draw ticks on the top side
                      #'ytick.right'         : True,   # draw ticks on the bottom side  
                      'xtick.major.size'    : 4,
                      'xtick.major.width'   : 1,
                      'xtick.minor.size'    : 2,
                      'xtick.minor.width'   : 1,
                      'ytick.major.size'    : 4,
                      'ytick.major.width'   : 1,
                      'ytick.minor.size'    : 2,
                      'ytick.minor.width'   : 1, 
                      'text.usetex'         : True,
                      'text.latex.unicode'  : True
                       }
            plt.rcParams.update(params)
            #plt.rc('xtick', labelsize=20)

            # Format axes
            nullfmt        = NullFormatter() 
            left, width    = 0.05, 0.9                                                                     #|These determine where the subplots go
            bottom, height = 0.1, 0.85
            #left_h = left+width+0.015
            #bottom_h = left+height+0.015

            box_centre    = [left, bottom, width, height]
                    
            # initialize figure
            fig = plt.figure(1, figsize=(12,8))
            plt.rc('xtick', labelsize=20)

            # Initialize subplots
            ax1 = plt.axes(box_centre)
            ax1.set_xlabel(r'Frequency [MHz]')
          
            bandwidth= (cfg_par['rfi']['highfreq'] - cfg_par['rfi']['lowfreq'] )/1e6

            if bandwidth <= 150.:
                step_bin=25.
            elif bandwidth <= 500. and bandwidth >  150.:
                step_bin=50.                
            else:
                step_bin=100.

            # Smaller multiple 
            freq_max = (freqs[-1] // step_bin) * step_bin  
            # Larger multiple 
            freq_min= (freqs[0] // step_bin) * step_bin + step_bin

            freqs_plot_bin=np.arange(freq_min,freq_max+step_bin,step_bin)
            freqs_plot_bin=np.round(freqs_plot_bin,0)
            freqs_plot_bin=np.array(freqs_plot_bin,dtype=int)

            #for i in xrange(0, len(freqs_plot_bin)):
            #    idx = (np.abs(freqs_plot_bin[i] - freqs)).argmin()
            #    freqs_plot_idx[i]=idx
            ax1.set_xticks(freqs_plot_bin)
            ax1.set_xticklabels(freqs_plot_bin)
            ax1.set_xlim([(cfg_par['rfi']['lowfreq']-5.*cfg_par['rfi']['chan_widths'])/1e6,
                (cfg_par['rfi']['highfreq']+5.*cfg_par['rfi']['chan_widths'])/1e6])
           
            if cfg_par['plots']['plot_noise'] != 'rfi':
                ax1.set_yscale('log', basey=10)
                if np.isnan(np.sum(noise_all)):
                    noise_all=np.zeros([len(freqs)])+100.
                if np.isnan(np.sum(noise_short)):
                    noise_short=np.zeros([len(freqs)])+100.
                if np.isnan(np.sum(noise_long)):
                    noise_long= np.zeros([len(freqs)])+100.

            #define title output                     
            
            #plot
            label_all = 'All baselines' 
            label_long = r'Baselines $>$ '+str(cfg_par['rfi']['baseline_cut'])+' m'
            label_short = r'Baselines $<$ '+str(cfg_par['rfi']['baseline_cut'])+' m' 

            if cfg_par['plots']['plot_long_short'] == True:
                self.logger.info("\t ... Plotting RFI in long and short baselines ...")
                ax1.step(freqs,noise_short, where= 'pre', color='red', linestyle='-',label=label_short)
                ax1.step(freqs,noise_long, where= 'pre', color='blue', linestyle='-',label=label_long)
                out_plot = out_plot+'_sl'

            ax1.step(freqs,noise_all, where= 'pre', color='black', linestyle='-',label=label_all)

            #titleplot = self.target+': '+self.aperfi_startime+' - '+self.aperfi_endtime
            #ax1.set_title(titleplot)
            
            # set axis, legend ticks

            #ax1.set_xlim([np.min(freqs)-5,np.max(freqs)+5])

            #xticks_num = np.linspace(cfg_par['rfi']['lowfreq'],cfg_par['rfi']['highfreq'],10,dtype=int)
            #ax1.set_xticks(xticks_num)

            if cfg_par['plots']['plot_noise']  == 'noise_factor':
                ax1.set_yticks([1,round(np.sqrt(2),2),2,3,5,10,50]) 
                ax1.set_yscale('linear')     
                ax1.set_ylabel(r'Factor of noise increase')
                ax1.get_yaxis().set_major_formatter(ticker.ScalarFormatter())

            if cfg_par['plots']['plot_noise'] == 'noise':
                ax1.set_yscale('linear')     
                ax1.set_ylabel(r'Predicted noise [mJy beam$^{-1}$]')     
                ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

            if cfg_par['plots']['plot_noise'] == 'rfi':
                if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                    ax1.set_ylabel(r'$\% > 5 \times$ r.m.s.')
                if cfg_par['rfi']['RFInder_mode'] == 'use_flags':
                    ax1.set_ylabel(r'$\%$ flagged visibilites')    
                ax1.set_ylim([-5,105]) 
                ax1.set_yticks([0,20,40,60,80,100]) 
            
            legend = plt.legend(loc=1)
            legend.get_frame().set_edgecolor('black')

            if cfg_par['rfi']['chunks']['time_enable']== False:
                start = cfg_par['rfi']['startdate']
                end = cfg_par['rfi']['enddate']

            elif cfg_par['rfi']['chunks']['time_enable'] == True:
                time_del = TimeDelta(time_delta*60., format='sec')
                time_delta_plus = TimeDelta(float(cfg_par['rfi']['chunks']['time_step'])*60., format='sec')
                start = cfg_par['rfi']['startdate']+time_del
                end = start+time_delta_plus
            
            pol =str(cfg_par['rfi']['polarization']) 
            if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                rfi_clip = str(cfg_par['rfi']['rms_clip'])+r'$\sigma$ clip'        
                title_plot = 'Stokes: {0:s} / {1:s} / {2:%d}{2:%b}{2:%y}: {2:%H}:{2:%M} - {3:%H}:{3:%M}'.format(pol,rfi_clip,start.datetime,end.datetime)
            
            if cfg_par['rfi']['RFInder_mode'] == 'use_flags':
                title_plot = 'Stokes: {0:s} / {1:s} / {2:%d}{2:%b}{2:%y}: {2:%H}:{2:%M} - {3:%H}:{3:%M}'.format(pol,'Flags',start.datetime,end.datetime)
            
            ax1.set_title(title_plot)
            ax1.minorticks_on()
            ax1.minorticks_on()
            # Save figure to file
            if cfg_par['rfi']['RFInder_mode']== 'use_flags':
                rfi_freq_plot = out_plot+'_flags.png'
            if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
                rfi_freq_plot = out_plot+'_rfi.png'

            plt.savefig(rfi_freq_plot,format='png',overwrite = True)      
            plt.close(fig)        
            self.logger.info("\t ... RFI in 1D plotted ...\n\n")
       

    def plot_altaz(self,cfg_par,number_chunks):
        '''
        Plots the elevation/azimuth of the observation scan binned by time chunks, for every binned spectral window
        '''


        self.logger.info("\t ... Plotting Alt/Az for binned dataset ... \n")
      
        if os.path.exists(cfg_par['general']['timetabledir']) == False:
            self.logger.error("\t Folder with time subsets missing")

        table_tmp = string.split(cfg_par['general']['msname'][0],'.MS')
        if len(table_tmp) == 1:
            table_tmp = string.split(cfg_par['general']['msname'][0],'.ms')

        freq_S = cfg_par['rfi']['lowfreq']/1e6
        freq_E = cfg_par['rfi']['highfreq']/1e6

        step_bin = cfg_par['rfi']['chunks']['spw_width']
        
        freqs_bin=np.arange(freq_S,freq_E+step_bin,step_bin)

        pol = cfg_par['rfi']['polarization'] 

        for j in xrange(0,len(freqs_bin)):
            spw = np.full((number_chunks),np.nan)
            azimuth = np.full((number_chunks),np.nan)
            altitude= np.full((number_chunks),np.nan)
            flags = np.full((number_chunks),np.nan)

            for i in xrange(0,number_chunks):
            
                tabledir = cfg_par['general']['timetabledir'] 
                time_delta = float(cfg_par['rfi']['chunks']['time_step'])*i
                time_tmp = int(time_delta)

                #name file according to WSRT beams
                if time_tmp == 0:
                    time_name = '00'+str(time_tmp)+'m'
                elif time_tmp <100:
                    time_name = '0'+str(time_tmp)+'m'
                else:
                    time_name= str(time_tmp)+'m'    

                if cfg_par['rfi']['RFInder_mode']== 'use_flags':
                    table_name = str(table_tmp[0])+'_flags_'+time_name+'_spwbin.fits'
                elif cfg_par['rfi']['RFInder_mode']== 'rms_clip':
                    table_name = str(table_tmp[0])+'_rfi_'+time_name+'_spwbin.fits'

                rfi_table = tabledir+table_name

                #open file
                if os.path.exists(rfi_table) == False:
                    self.logger.error('### Table of RFI results does not exist ###')    
                    continue         
                else:    
                    table = Table.read(rfi_table)
                    spw[i] = table['frequency'][j]
                    azimuth[i] = table['azimuth'][j]
                    altitude[i] = table['altitude'][j]
                    flags[i] = table['percentage_flags'][j]
            #az=np.array([az],dtype=float)
            #alt=np.array([alt],dtype=float)
            #az,alt=self.aitoff(az,alt)
            #azimuth=np.array(azimuth,dtype=float)
            #altitude=np.array(altitude,dtype=float)

            azimuth,altitude =self.aitoff(azimuth-180.,altitude)
            #azimuth=np.array([azimuth],dtype=float)
            #altitude=np.array([altitude],dtype=float)

            plotdir = cfg_par['general']['altazplotdir']

            start_freq = int(np.round(freq_S+float(cfg_par['rfi']['chunks']['spw_width'])*j,0))
            end_freq = int(np.round(start_freq+float(cfg_par['rfi']['chunks']['spw_width']),0))

            spwname= str(start_freq)+'-'+str(end_freq)

            # Save figure to file
            if cfg_par['rfi']['RFInder_mode']== 'use_flags':
               altazplot = plotdir+'AltAZ_flags'+spwname+'MHz.png'
            if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
               altazplot = plotdir+'AltAZ_rfi'+spwname+'MHz.png'

                    # initialize plotting parameters
            params = {'font.family'         :' serif',
                      'font.style'          : 'normal',
                      'font.weight'         : 'book',
                      'font.size'           : 18.0,
                      'axes.linewidth'      : 1,
                      'lines.linewidth'     : 1,
                      'xtick.labelsize'     : 16,
                      'ytick.labelsize'     : 16, 
                      'xtick.direction'     :'in',
                      'ytick.direction'     :'in',
                      #'xtick.top'           : True,   # draw ticks on the top side
                      #'xtick.bottom'        : True,   # draw ticks on the bottom side    
                      #'ytick.left'          : True,   # draw ticks on the top side
                      #'ytick.right'         : True,   # draw ticks on the bottom side  
                      'xtick.major.size'    : 4,
                      'xtick.major.width'   : 1,
                      'xtick.minor.size'    : 2,
                      'xtick.minor.width'   : 1,
                      'ytick.major.size'    : 4,
                      'ytick.major.width'   : 1,
                      'ytick.minor.size'    : 2,
                      'ytick.minor.width'   : 1, 
                      'text.usetex'         : True,
                      'text.latex.unicode'  : True
                       }
            plt.rcParams.update(params)
            #plt.rc('xtick', labelsize=20)
            bandwidth= (cfg_par['rfi']['highfreq'] - cfg_par['rfi']['lowfreq'] )/1e6

            if bandwidth <= 500.:
                colormap = 'nipy_spectral_r'
            else:
                colormap = 'afmhot'

            # Format axes
            nullfmt        = NullFormatter() 
            left, width    = 0.12, 0.64                                                                     #|These determine where the subplots go
            bottom, height = 0.12, 0.54
            left_h = left+width+0.015
            bottom_h = left+height+0.015

            box_centre    = [left, bottom, width, height]
            box_x      = [left, bottom_h, width, 0.12]
            box_y      = [left_h, bottom, 0.12, height]
            box_cbar       = [left_h+0.13, bottom, 0.05, height]

            fig = plt.figure(1, figsize=(12,8))

            ax_centre = plt.axes(box_centre)
            ax_x = plt.axes(box_x)
            ax_y = plt.axes(box_y)
            ax_cbar   = plt.axes(box_cbar)

            # Remove labels from supplementary plots
            ax_x.xaxis.set_major_formatter(nullfmt)
            ax_y.yaxis.set_major_formatter(nullfmt)
            ax_cbar.xaxis.set_major_formatter(nullfmt)
            ax_cbar.yaxis.set_major_formatter(nullfmt)
            ax_cbar.set_visible(False)

            # Add minor tick marks
            ax_centre.minorticks_on()
            ax_x.minorticks_on()
            ax_y.minorticks_on()

            # initialize figure
            #fig = plt.figure(figsize =(8,8))
            #fig.subplots_adjust(hspace=0.0)
            #gs = gridspec.GridSpec(1, 1)

            # Initialize subplots
            #ax1 = fig.add_subplot(gs[0])
            #centre
            ax_centre.set_xlabel(r'Azimuth [deg]')
            ax_centre.set_ylabel(r'Altitude [deg]',labelpad=15)
            
            ax_centre.set_ylim([0,90])
            ax_centre.set_xlim([-180.,180])
            ax_centre.set_xticks([-180,-135,-90,-45,0,45,90,135,180])        
            #ax_centre.set_xticklabels(['0','45','90','135','180','225','270','315','0'])
            #ax_centre.set_yticks([0,10,20,30,40,50,60,70,80,90])
            #ax_centre.set_xticks([])        
            ax_centre.set_yticks([])

            asse = ax_centre.scatter(azimuth,altitude,c=flags,cmap=colormap,vmin=0,vmax=100.,s=100,linewidth=0.5,edgecolors='black',label=r'Stokes: {0:s}'.format(pol))

            azs=np.arange(-180.,181.,45.)
            als=np.arange(0.,91.,1.)
            azs[azs==0]=1e-6
            als[als==0]=1e-6
            for counter in azs:
                [X,Y]=self.aitoff(counter*np.ones(als.shape),als)
                ax_centre.plot(X,Y,'k-',linewidth=0.5)

            azs=np.arange(-180.,181.,1.)
            als=np.arange(0.,91.,15.)
            azs[azs==0]=1e-6
            als[als==0]=1e-6
            for al in als:
                [X,Y]=self.aitoff(azs,al*np.ones(azs.shape))
                ax_centre.plot(X,Y,'k-',linewidth=0.5)

            for dd in [15.,30.,45.,60.,75.]:
                [X,Y]=self.aitoff(-180.,dd)
                ax_centre.text(X,Y,r'${0:s}$'.format(str(int(round(dd,0)))),horizontalalignment='right',fontsize='medium')

            legend = ax_centre.legend(loc=1,fontsize=14,handletextpad=0.1,borderpad=0.2)
            legend.get_frame().set_edgecolor('black')
            legend.legendHandles[0].set_color('black')
            legend.legendHandles[0]._sizes = [35]

            start = cfg_par['rfi']['startdate']
            end = cfg_par['rfi']['enddate']
            
             # colorbar
            cbar = plt.colorbar(asse,ax=ax_cbar, fraction=1.0,ticks=[0,20,40,60,80,100]) 
            if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                cbar.set_label(r'$\% > 5 \times$ r.m.s.')
            if cfg_par['rfi']['RFInder_mode'] == 'use_flags':
                cbar.set_label(r'$\%$ flagged visibilites')       
            #title
            if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                rfi_clip = str(cfg_par['rfi']['rms_clip'])+r'$\sigma$ clip'        
                title_plot = '{0:d}-{1:d} MHz / {2:s} / {3:%d}{3:%b}{3:%y}: {3:%H}:{3:%M} - {4:%H}:{4:%M}'.format(start_freq,end_freq,rfi_clip,start.datetime,end.datetime)
            if cfg_par['rfi']['RFInder_mode'] == 'use_flags':
                title_plot = '{0:d}-{1:d} MHz / {2:s} / {3:%d}{3:%b}{3:%y}: {3:%H}:{3:%M} - {4:%H}:{4:%M}'.format(start_freq,end_freq,'Flags',start.datetime,end.datetime)        

            #x plot
            #ax_x.set_xlabel(r'Azimuth [deg]',fontsize=16)
            if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                ax_x.set_ylabel(r'$\%$ RFI')
            if cfg_par['rfi']['RFInder_mode'] == 'use_flags':
                ax_x.set_ylabel(r'$\%$ flagged')                
            ax_x.set_ylim([-5,100])
            ax_x.set_xlim([-180,180])
            ax_x.set_yticks([0,25,50,75,100])
   
            ax_x.set_title(title_plot)
           
            ax_x.scatter(azimuth,flags,c=flags,cmap=colormap,vmin=0,vmax=100.,s=70,linewidth=0.5,edgecolors='black')
            
            #y plot
            #ax_x.set_xlabel(r'Azimuth [deg]',fontsize=16)
            if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
                ax_y.set_xlabel(r'$\%$ RFI ')
            if cfg_par['rfi']['RFInder_mode'] == 'use_flags':
                ax_y.set_xlabel(r'$\%$ flagged')          

            ax_y.set_xlim([-5,100])
            ax_y.set_ylim([0,90])
            ax_y.set_xticks([0,25,50,75,100])
            ax_y.scatter(flags,altitude,c=flags,cmap=colormap,vmin=0,vmax=100.,s=70,linewidth=0.5,edgecolors='black')


            # Finish everything up
            plt.savefig(altazplot,format='png',overwrite = True)
            plt.close()
            self.logger.info(('\t ... ALT/AZ for spw: {0:d}-{1:d} MHz  ...\n').format(start_freq,end_freq))

        return 0

    def plot_altaz_short(self,cfg_par):
        '''
        Plots the elevation/azimuth of the observation scan binned by time chunks, for every binned spectral window
        '''


        self.logger.info("\t ... Plotting Alt/Az for observation ... \n")
        
        times = cfg_par['rfi']['times']
        azimuth = []
        altitude = []

        for i in xrange(0,len(times)):
        
            altaz = rfiST.alt_az(cfg_par,times[i])
            azimuth.append(altaz.az)
            altitude.append(altaz.alt)


        azimuth=np.array([azimuth],dtype=float)
        altitude=np.array([altitude],dtype=float)

        [azimuth,altitude] =self.aitoff(azimuth-180.,altitude)


        plotdir = cfg_par['general']['plotdir']

        freq_S = cfg_par['rfi']['lowfreq']/1e6
        freq_E = cfg_par['rfi']['highfreq']/1e6

        start_freq = int(np.round(freq_S,0))
        end_freq = int(np.round(freq_E,0))

        spwname= str(start_freq)+'-'+str(end_freq)

        # Save figure to file
        altazplot = cfg_par['general']['plotdir']+'AltAZ_full.png'

                # initialize plotting parameters
        params = {'font.family'         :' serif',
                  'font.style'          : 'normal',
                  'font.weight'         : 'book',
                  'font.size'           : 18.0,
                  'axes.linewidth'      : 1,
                  'lines.linewidth'     : 1,
                  'xtick.labelsize'     : 16,
                  'ytick.labelsize'     : 16, 
                  'xtick.direction'     :'in',
                  'ytick.direction'     :'in',
                  'xtick.major.size'    : 4,
                  'xtick.major.width'   : 1,
                  'xtick.minor.size'    : 2,
                  'xtick.minor.width'   : 1,
                  'ytick.major.size'    : 4,
                  'ytick.major.width'   : 1,
                  'ytick.minor.size'    : 2,
                  'ytick.minor.width'   : 1, 
                  'text.usetex'         : True,
                  'text.latex.unicode'  : True
                   }
        plt.rcParams.update(params)

        # Format axes
        nullfmt        = NullFormatter() 
        left, width    = 0.05, 0.9                                                                     #|These determine where the subplots go
        bottom, height = 0.1, 0.85

        box_centre    = [left, bottom, width, height]
        fig = plt.figure(1, figsize=(12,8))

        ax_centre = plt.axes(box_centre)
        # Add minor tick marks
        ax_centre.minorticks_on()

        ax_centre.set_xlabel(r'Azimuth [deg]')
        ax_centre.set_ylabel(r'Altitude [deg]')
        
        ax_centre.set_ylim([0,90])
        ax_centre.set_xlim([-180.,180])
        ax_centre.set_xticks([-180,-135,-90,-45,0,45,90,135,180])        
        ax_centre.set_yticks([])

        asse = ax_centre.scatter(azimuth,altitude,c='black',s=100,linewidth=0.5,edgecolors='black',label=r'${0:s}$'.format(cfg_par['general']['fieldname']))

        azs=np.arange(-180.,181.,45.)
        als=np.arange(0.,91.,1.)
        azs[azs==0]=1e-6
        als[als==0]=1e-6
        for counter in azs:
            [X,Y]=self.aitoff(counter*np.ones(als.shape),als)
            ax_centre.plot(X,Y,'k-',linewidth=0.5)

        azs=np.arange(-180.,181.,1.)
        als=np.arange(0.,91.,15.)
        azs[azs==0]=1e-6
        als[als==0]=1e-6
        for al in als:
            [X,Y]=self.aitoff(azs,al*np.ones(azs.shape))
            ax_centre.plot(X,Y,'k-',linewidth=0.5)

        for dd in [15.,30.,45.,60.,75.]:
            [X,Y]=self.aitoff(-180.,dd)
            ax_centre.text(X,Y,r'${0:s}$'.format(str(int(round(dd,0)))),horizontalalignment='right',fontsize='medium')

        legend = ax_centre.legend(loc=1,fontsize=14,handletextpad=0.1,borderpad=0.2)
        legend.get_frame().set_edgecolor('black')
        legend.legendHandles[0].set_color('black')
        legend.legendHandles[0]._sizes = [35]

        start = cfg_par['rfi']['startdate']
        end = cfg_par['rfi']['enddate']

        title_plot = '{0:d}-{1:d} MHz / {2:%d}{2:%b}{2:%y}: {2:%H}:{2:%M} - {3:%H}:{3:%M}'.format(start_freq,end_freq,start.datetime,end.datetime)
         
        ax_centre.set_title(title_plot)

        # Finish everything up
        plt.savefig(altazplot,format='png',overwrite = True)
        plt.close()
        self.logger.info(('\t ... ALT/AZ for {0:s} ...\n').format(cfg_par['general']['fieldname']))

        return 0
    

    def gif_me_up(self,cfg_par,filenames,outmovie):
        

        self.logger.info(('\t ... Creating movie ...'))
       
        # initiate an empty  list of "plotted" images 
        images=[]

        for filename in filenames:
            images.append(io.imread(filename))
            
        io.mimsave(outmovie,images,duration=1.5)


        self.logger.info(('\t ... Movie written to file ...\n'))

        return 0


