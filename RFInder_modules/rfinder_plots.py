import os,string,sys
import matplotlib
#matplotlib.use('Qt5Agg')
#matplotlib.rcParams['backend']='template'
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import ticker

import aplpy
from astropy.io import fits as fits
import numpy as np
import os
from astropy.table import Table
import pyrap.tables as tables
from astropy.time import Time, TimeDelta
import rfi
rfi = rfi.rfi()

import logging


logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('PLOT')


def plot_rfi_im(cfg_par,time_step=-1):      
    '''
    
    Plots the .fits image output of rfi_im jpg format.
    
    '''        
    
    #check if image exists
    #plot image
    logger.info("\t ... Plotting RFI in 2D ... \n")

    if time_step != -1:
        outputdir = str(cfg_par['general']['rfitimedir'])
        plotdir = cfg_par['general']['timeplotdir']
        time_tmp = int(float(cfg_par['rfi']['chunks']['time_step'])*time_step)
        if time_tmp == 0:
            time_name = '00'+str(time_tmp)+'m'
        elif time_tmp <100:
            time_name = '0'+str(time_tmp)+'m'
        else:
            time_name= str(time_tmp)+'m'    
    else:
        outputdir = str(cfg_par['general']['rfidir'])
        plotdir = str(cfg_par['general']['plotdir'])        
        time_name = 'full' 

    rfi_freq_base = outputdir+'freq_base_'+time_name+'_im.fits'   
    rfi_freq_base_plot = plotdir+'freq_base_'+time_name+'_pl.png'
    


    fig = aplpy.FITSFigure(rfi_freq_base,figsize=(12,8))

    #plot colorscale & colorbar
    fig.show_colorscale(aspect='auto', cmap='nipy_spectral_r',vmin=0,vmax=100)
    fig.show_colorbar()
    fig.colorbar.set_width(0.2)
    fig.colorbar.set_font(size=20, weight='medium', \
                          stretch='normal', family='sans-serif', \
                          style='normal', variant='normal')
    fig.colorbar.set_axis_label_font(size=20)
    fig.colorbar.set_axis_label_text(r'$\% > 5 \times$ r.m.s.')

    #set axis
    fig.axis_labels.set_font(size=20, weight='medium', \
                             stretch='normal', family='sans-serif', \
                             style='normal', variant='normal')
    fig.tick_labels.set_font(size=20, weight='medium', \
                             stretch='normal', family='sans-serif', \
                             style='normal', variant='normal') 
    #titleplot = self.target+': '+self.aperfi_startime+' - '+self.aperfi_endtime

    plt.savefig(rfi_freq_base_plot,format='png' ,overwrite=True)

    logger.info("\t ... RFI per baseline lenght and frequency plotted ... \n\n")

def plot_rfi_imshow(cfg_par,time_step=-1):      
    '''
    
    Plots the .fits image output of rfi_im jpg format.
    
    '''        
    
    #check if image exists
    #plot image

    logger.info("\t ... Plotting RFI in 2D ... \n")


    if time_step != -1:
        outputdir = cfg_par['general']['rfitimedir']
        plotdir = cfg_par['general']['timeplotdir']
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

    if cfg_par['rfi']['use_flags']== True:
        rfi_freq_base =outputdir+'freq_base_'+time_name+'_flags.fits'
        rfi_freq_base_plot = plotdir+'freq_base_'+time_name+'_flags.png'
    if cfg_par['rfi']['use_flags']== False:
        rfi_freq_base =outputdir+'freq_base_'+time_name+'_rfi.fits'
        rfi_freq_base_plot = plotdir+'freq_base_'+time_name+'_rfi.png'
    

    # open file
    t = fits.open(rfi_freq_base)
    data = t[0].data
    prihdr= t[0].header
    freqs = (np.linspace(1, data.shape[1], data.shape[1])- prihdr['CRPIX1'])*prihdr['CDELT1'] + prihdr['CRVAL1']
    
    # set x-y axes
    step_bin=50.
    freqs_plot_bin=np.arange(freqs[0],freqs[-1]+step_bin,step_bin)
    freqs_plot_bin=np.round(freqs_plot_bin,0)

    freqs_plot_idx=np.zeros(freqs_plot_bin.shape)
    for i in xrange(0, len(freqs_plot_bin)):
        idx = (np.abs(freqs_plot_bin[i] - freqs)).argmin()
        freqs_plot_idx[i]=idx


    tele= cfg_par['rfi']['telescope']
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
    
    # initialize plotting parameters
    params = {'font.family'         :' serif',
              'font.style'          : 'normal',
              'font.weight'         : 'medium',
              'font.size'           : 20.0,
              'text.usetex': True,
              'text.latex.unicode': True
               }
    plt.rcParams.update(params)


    fig, ax = plt.subplots(figsize=(12,8))
    im = ax.imshow(data,vmin=0,vmax=100,cmap='nipy_spectral_r')

    # ticks & labels
    ax.set_xticks(freqs_plot_idx)
    ax.set_yticks(input_baselines_idx)

    ax.set_xticklabels(freqs_plot_bin)
    ax.set_yticklabels(input_baselines)

    ax.set_xlabel('Frequency [MHz]')
    ax.set_ylabel('Baseline lenght [m]')
    
    # colorbar
    cbar = fig.colorbar(im,ticks=[10, 20, 30,40,50,60,70,80,90,100]) 
    if cfg_par['rfi']['use_flags'] == False:
        cbar.set_label(r'$\% > 5 \times$ r.m.s.',size=20)
    if cfg_par['rfi']['use_flags'] == True:
        cbar.set_label(r'$\%$ flagged visibilites',size=20)

    #times, start, end = rfi.time_chunk(cfg_par)
    #start.format='iso' 
    #start.subformat='date_hm'   
    #end.format='iso'
    #end.subformat='date_hm'


    if cfg_par['rfi']['chunks']['time_enable']== False:
        start = cfg_par['rfi']['startdate']
        end = cfg_par['rfi']['enddate']

    elif cfg_par['rfi']['chunks']['time_enable'] == True:
        time_delta = dt2 = TimeDelta(time_delta*60., format='sec')
        start = cfg_par['rfi']['startdate']
        end = start+ time_delta

    if cfg_par['rfi']['use_flags'] == False:
        title_plot = '{0:s} / {1:%d}{1:%b}{1:%y}: {1:%H}:{1:%M} - {2:%H}:{2:%M}'.format('RFI clip',start.datetime,end.datetime)
    if cfg_par['rfi']['use_flags'] == True:
        title_plot = '{0:s} / {1:%d}{1:%b}{1:%y}: {1:%H}:{1:%M} - {2:%H}:{2:%M}'.format('Flags',start.datetime,end.datetime)
    
    ax.set_title(title_plot)
 
    plt.savefig(rfi_freq_base_plot,format='png' ,overwrite=True)

    logger.info("\t ... RFI per baseline lenght and frequency plotted ... \n\n")

def plot_noise_frequency(cfg_par,time_step=-1):
    '''
    Plots the noise or % of rfi per frequency channel for all, long and short baselines.
    In default.cfga
        aperfi_noise = 'rfi' (or flag)
        aperfi_plot_long_short = False
    '''

    #open file
    #if os.path.exists(self.rfi_table) == False:
    #    self.logger.error('### Table of RFI and flags of visibilities does not exist ###')    
    #    self.logger.error('### Run aperfi.rfi_frequency() first ###')  
    #else:  

    logger.info("\t ... Plotting RFI in 1D ... \n")


    table_tmp = string.split(cfg_par['general']['msname'][0],'.MS')
    if len(table_tmp) == 1:
        table_tmp = string.split(cfg_par['general']['msname'][0],'.ms')

    if time_step != -1:
        tabledir = cfg_par['general']['timetabledir'] 
        plotdir = cfg_par['general']['timeplotdir']
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

    if cfg_par['rfi']['use_flags']== True:
        table_name = str(table_tmp[0])+'_'+time_name+'_flags.fits'
    elif cfg_par['rfi']['use_flags']== False:
        table_name = str(table_tmp[0])+'_'+time_name+'_rfi.fits'

    rfi_table = tabledir+table_name
    

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
        rms = np.array(cfg_par['rfi']['theo_rms'],dtype=float)
        noise_all = noise_factor*rms
        noise_short = noise_factor_short*rms
        noise_long = noise_factor_long*rms
        out_plot = plotdir+'freq_noise_'+time_name
    
    if cfg_par['plots']['plot_noise'] == 'noise_factor':
        logger.info("\t ... Plotting factor of noise increas per frequency channel ...")
        noise_all = noise_factor
        noise_short = noise_factor_short
        noise_long = noise_factor_long    
        out_plot = plotdir+'freq_noise_factor_'+time_name
    
    if cfg_par['plots']['plot_noise'] == 'rfi':
        logger.info("\t ... Plotting percentage of flagged RFI per frequency channel ...")
        noise_all = flags
        noise_long = flags_long
        noise_short = flags_short
        out_plot = plotdir+'freq_flags_'+time_name


    # initialize plotting parameters
    params = {'font.family'         :' serif',
              'font.style'          : 'normal',
              'font.weight'         : 'medium',
              'font.size'           : 20.0,
              'text.usetex': True,
              'text.latex.unicode': True
               }
    plt.rcParams.update(params)
    
    # initialize figure
    fig = plt.figure(figsize =(14,8))
    fig.subplots_adjust(hspace=0.0)
    gs = gridspec.GridSpec(1, 1)
    plt.rc('xtick', labelsize=20)

    # Initialize subplots
    ax1 = fig.add_subplot(gs[0])
    ax1.set_xlabel(r'Frequency [MHz]',fontsize=20)
    
    if cfg_par['plots']['plot_noise'] != 'rfi':
        ax1.set_yscale('log', basey=10)

    #define title output                     
    
    #plot
    label_all = 'All baselines' 
    label_long = r'Baselines $>$ '+str(cfg_par['rfi']['baseline_cut'])+' m'
    label_short = r'Baselines $<$ '+str(cfg_par['rfi']['baseline_cut'])+' m' 

    if cfg_par['plots']['plot_long_short'] == True:
        logger.info("\t ... Plotting RFI in long and short baselines ...")
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
        if cfg_par['rfi']['use_flags'] == False:
            ax1.set_ylabel(r'$\% > 5 \times$ r.m.s.',size=20)
        if cfg_par['rfi']['use_flags'] == True:
            ax1.set_ylabel(r'$\%$ flagged visibilites',size=20)    
    
    legend = plt.legend()
    legend.get_frame().set_edgecolor('black')

    if cfg_par['rfi']['chunks']['time_enable']== False:
        start = cfg_par['rfi']['startdate']
        end = cfg_par['rfi']['enddate']

    elif cfg_par['rfi']['chunks']['time_enable'] == True:
        time_delta = dt2 = TimeDelta(time_delta*60., format='sec')
        start = cfg_par['rfi']['startdate']
        end = start+ time_delta

    if cfg_par['rfi']['use_flags'] == False:
        title_plot = '{0:s} / {1:%d}{1:%b}{1:%y}: {1:%H}:{1:%M} - {2:%H}:{2:%M}'.format('RFI clip',start.datetime,end.datetime)
    if cfg_par['rfi']['use_flags'] == True:
        title_plot = '{0:s} / {1:%d}{1:%b}{1:%y}: {1:%H}:{1:%M} - {2:%H}:{2:%M}'.format('Flags',start.datetime,end.datetime)
    
    ax1.set_title(title_plot)

    # Save figure to file
    if cfg_par['rfi']['use_flags']== True:
        rfi_freq_plot = out_plot+'_flags.png'
    if cfg_par['rfi']['use_flags']== False:
        rfi_freq_plot = out_plot+'_rfi.png'

    plt.savefig(rfi_freq_plot,format='png',overwrite = True)      

    logger.info("\t ... RFI in 1D plotted ...\n\n")
   

def plot_altaz(cfg_par,number_chunks):
    '''
    Plots the elevation/azimuth of the observation scan binned by time chunks, for every binned spectral window
    '''


    logger.info("\t ... Plotting Alt/Az for binned dataset ... \n")
  


    if os.path.exists(cfg_par['general']['timetabledir']) == False:
        logger.error("\t Folder with time subsets missing")

    table_tmp = string.split(cfg_par['general']['msname'][0],'.MS')
    if len(table_tmp) == 1:
        table_tmp = string.split(cfg_par['general']['msname'][0],'.ms')

    freq_S = cfg_par['rfi']['lowfreq']/1e6
    freq_E = cfg_par['rfi']['highfreq']/1e6

    step_bin = cfg_par['rfi']['chunks']['spw_width']
    
    freqs_bin=np.arange(freq_S,freq_E+step_bin,step_bin)


    for j in xrange(0,len(freqs_bin)):
        spw = []
        az = []
        alt = []
        flags =[]


        for i in xrange(0,number_chunks):
        
            tabledir = cfg_par['general']['timetabledir'] 
            plotdir = cfg_par['general']['plotdir']
            time_delta = float(cfg_par['rfi']['chunks']['time_step'])*i
            time_tmp = int(time_delta)

            #name file according to WSRT beams
            if time_tmp == 0:
                time_name = '00'+str(time_tmp)+'m'
            elif time_tmp <100:
                time_name = '0'+str(time_tmp)+'m'
            else:
                time_name= str(time_tmp)+'m'    

            if cfg_par['rfi']['use_flags']== True:
                table_name = str(table_tmp[0])+'_'+time_name+'_spwbin_flags.fits'
            elif cfg_par['rfi']['use_flags']== False:
                table_name = str(table_tmp[0])+'_'+time_name+'_spwbin_rfi.fits'

            rfi_table = tabledir+table_name

            table = Table.read(rfi_table)

            spw.append(table['frequency'][j])
            az.append(table['azimuth'][j])
            alt.append(table['elevation'][j])
            flags.append(table['percentage_flags'][j])

        plotdir = cfg_par['general']['altazplotdir']

        start_freq = int(np.round(freq_S+float(cfg_par['rfi']['chunks']['spw_width'])*j,0))
        end_freq = int(np.round(start_freq+float(cfg_par['rfi']['chunks']['spw_width']),0))

        spwname= str(start_freq)+'-'+str(end_freq)
        # Save figure to file
        if cfg_par['rfi']['use_flags']== True:
           altazplot = plotdir+'AltAZ_'+spwname+'MHZ_flags.png'
        if cfg_par['rfi']['use_flags']== False:
           altazplot = plotdir+'AltAZ_'+spwname+'_rfi.png'


        start = cfg_par['rfi']['startdate']
        end = cfg_par['rfi']['enddate']


        # initialize plotting parameters
        params = {'font.family'         :' serif',
                  'font.style'          : 'normal',
                  'font.weight'         : 'medium',
                  'font.size'           : 18.0,
                  'text.usetex': True,
                  'text.latex.unicode': True
                   }
        plt.rcParams.update(params)
        
        # initialize figure
        fig = plt.figure(figsize =(8,8))
        fig.subplots_adjust(hspace=0.0)
        gs = gridspec.GridSpec(1, 1)
        plt.rc('xtick', labelsize=20)

        # Initialize subplots
        ax1 = fig.add_subplot(gs[0])
        ax1.set_xlabel(r'Azimuth [$^{\circ}$]',fontsize=20)
        ax1.set_ylabel(r'Elevation [$^{\circ}$]',fontsize=20)
        ax1.set_ylim([0,90])
        ax1.set_xlim([0,360])
        
        asse = ax1.scatter(az,alt,c=flags,cmap='nipy_spectral_r',vmin=0,vmax=100.)

         # colorbar
        cbar = plt.colorbar(asse,ticks=[10, 20, 30,40,50,60,70,80,90,100]) 
        if cfg_par['rfi']['use_flags'] == False:
            cbar.set_label(r'$\% > 5 \times$ r.m.s.',size=20)
        if cfg_par['rfi']['use_flags'] == True:
            cbar.set_label(r'$\%$ flagged visibilites',size=20)       

        if cfg_par['rfi']['use_flags'] == False:
            title_plot = 'SPW {0:d}-{1:d} MHz / {2:s} / {3:%d}{3:%b}{3:%y}: {3:%H}:{3:%M} - {4:%H}:{4:%M}'.format(start_freq,end_freq,'RFI clip',start.datetime,end.datetime)
        if cfg_par['rfi']['use_flags'] == True:
            title_plot = 'SPW {0:d}-{1:d} MHz / {2:s} / {3:%d}{3:%b}{3:%y}: {3:%H}:{3:%M} - {4:%H}:{4:%M}'.format(start_freq,end_freq,'Flags',start.datetime,end.datetime)        

        ax1.set_title(title_plot)

        plt.savefig(altazplot,format='png',overwrite = True)      

        logger.info(('\t ... ALT/AZ for spw: {0:d}-{1:d} MHz  ...').format(start_freq,end_freq))
   

    return 0