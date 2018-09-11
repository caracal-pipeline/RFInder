import os,string,sys
import matplotlib
#matplotlib.use('Qt5Agg')
#matplotlib.rcParams['backend']='template'
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib import ticker

import aplpy
from astropy.io import fits as pyfits
import numpy as np
import os
import pyrap.tables as tables


import rfi
rfi = rfi.rfi()

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def rfi_frequency(cfg_par,time_step=-1):
    '''
    Determines the rfi per frequency channel. Saves results in table rfi_table.fits
    For each channel the flag % and factor of noise increase are stored for all, long and short baselines
    Long and short baselines are separated in half, depending on the number of baselines
    '''
    table_tmp = string.split(cfg_par['general']['msname'][0],'.MS')

    if time_step != -1:
        time_tmp = int(float(cfg_par['rfi']['time_chunks']['time_step'])*time_step)
        if time_tmp == 0:
            time_name = '00'+str(time_tmp)+'m'
        elif time_tmp <100:
            time_name = '0'+str(time_tmp)+'m'
        else:
            time_name= str(time_tmp)+'m'
    else:
        time_name = 'full'
    

    table_name = str(table_tmp[0])+'_'+time_name+'.fits'
    
    rfi_table = cfg_par['general']['tabledir']+table_name
    rfi_freq_base = cfg_par['general']['rfidir']+'freq_base_'+time_name+'_im.fits'   

    #open file
    #if os.path.exists(self.rfi_freq_base) == False:
    #    self.logger.error('### Image of RFI sorted by frequency over baseline lenght does not exist ###')    
    #    self.logger.error('### Run aperfi.rfi_im() first ###')  
    #else:    
        
    # read data and header
    hdulist = pyfits.open(rfi_freq_base)  # read input                
    datacube = hdulist[0].data    
    prihdr = hdulist[0].header

    #set array of frequencies
    freqs = (np.linspace(1, datacube.shape[1], datacube.shape[1])\
                 - prihdr['CRPIX1'])*prihdr['CDELT1'] + prihdr['CRVAL1']
    
    # set y-array
    rms_lin = np.zeros([datacube.shape[1]])    
    flag_lin = np.zeros([datacube.shape[1]])    
    rms_lin_long = np.zeros([datacube.shape[1]]) + np.sqrt(2.)          
    rms_lin_short = np.zeros([datacube.shape[1]]) + np.sqrt(2.)   
    flag_lin_long = np.zeros([datacube.shape[1]]) + 50.          
    flag_lin_short = np.zeros([datacube.shape[1]]) + 50.

    for i in xrange(0,datacube.shape[1]):
        
        flag_lin_tmp = np.divide(np.sum(datacube[:,i]),datacube.shape[0])
        flag_lin[i] = flag_lin_tmp


        baseline_cutoff = float(cfg_par['rfi']['baseline_cut'])
        lenghts = np.array([cfg_par['rfi']['baseline_lenghts']])+0.
        idx = (np.abs(lenghts - baseline_cutoff)).argmin()

        shortbase=datacube[:idx,i]
        longbase = datacube[idx:,i]               
        
        rms_lin_tmp = 1.-np.divide(np.divide(np.sum(datacube[:,i]),datacube.shape[0]),100.)
        rms_lin[i] = np.divide(1.,np.sqrt(rms_lin_tmp))

        flag_lin_tmp = np.divide(np.sum(shortbase),len(shortbase))
        flag_lin_short[i] = flag_lin_tmp
        rms_lin_tmp_short = 1.-np.divide(np.divide(np.sum(shortbase),len(shortbase)),100.)
        rms_lin_short[i] *= np.divide(1.,np.sqrt(rms_lin_tmp_short))

        flag_lin_tmp = np.divide(np.sum(longbase),len(longbase))
        flag_lin_long[i] = flag_lin_tmp
        rms_lin_tmp_long = 1.-np.divide(np.divide(np.sum(longbase),len(longbase)),100.)
        rms_lin_long[i] *= np.divide(1.,np.sqrt(rms_lin_tmp_long))

    
    #rebin results
    if cfg_par['rfi']['spw_average']['enable'] == True:

        table_name_bin = str(table_tmp[0])+'_'+time_name+'_spwbin.fits'
        rfi_table_bin = cfg_par['general']['tabledir']+table_name_bin

        step_bin = cfg_par['rfi']['spw_average']['spw_width']


        freqs_bin=np.arange(freqs[0],freqs[-1]+step_bin,step_bin)
        flag_lin_bin=np.zeros(freqs_bin.shape)
        rms_lin_bin=np.zeros(freqs_bin.shape)
        flag_lin_bin_short=np.zeros(freqs_bin.shape)
        rms_lin_bin_short=np.zeros(freqs_bin.shape)
        flag_lin_bin_long=np.zeros(freqs_bin.shape)
        rms_lin_bin_long=np.zeros(freqs_bin.shape)



        for i in xrange(0, len(freqs_bin)-1):
            #look for the right velocity bin
            index = (freqs_bin[i] <= freqs) & (freqs < freqs_bin[i+1])

            flag_lin_bin[i] = np.nanmean(flag_lin[index])
            rms_lin_bin[i] = np.nanmean(rms_lin[index])
            flag_lin_bin_short[i] = np.nanmean(flag_lin_short[index])
            rms_lin_bin_short[i] = np.nanmean(rms_lin_short[index])
            flag_lin_bin_long[i] = np.nanmean(flag_lin_long[index])
            rms_lin_bin_long[i] = np.nanmean(rms_lin_long[index])

        # save fits table        
        c1 = pyfits.Column(name='frequency', format='D', unit='MHz', array=freqs_bin)
        c2 = pyfits.Column(name='percentage_flags', format='D', unit='-', array=flag_lin_bin)
        c3 = pyfits.Column(name='noise_factor', format='D', unit = '-', array=rms_lin_bin)
        c4 = pyfits.Column(name='percentage_flags_short', format='D', unit='-', array=flag_lin_bin_short)
        c5 = pyfits.Column(name='noise_factor_short', format='D', unit = '-', array=rms_lin_bin_short)
        c6 = pyfits.Column(name='percentage_flags_long', format='D', unit='-', array=flag_lin_bin_long)
        c7 = pyfits.Column(name='noise_factor_long', format='D', array=rms_lin_bin_long)        

        fits_table = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7])    

        fits_table.writeto(rfi_table_bin, overwrite = True)


    # save fits table        
    c1 = pyfits.Column(name='frequency', format='D', unit='MHz', array=freqs)
    c2 = pyfits.Column(name='percentage_flags', format='D', unit='-', array=flag_lin)
    c3 = pyfits.Column(name='noise_factor', format='D', unit = '-', array=rms_lin)
    c4 = pyfits.Column(name='percentage_flags_short', format='D', unit='-', array=flag_lin_short)
    c5 = pyfits.Column(name='noise_factor_short', format='D', unit = '-', array=rms_lin_short)
    c6 = pyfits.Column(name='percentage_flags_long', format='D', unit='-', array=flag_lin_long)
    c7 = pyfits.Column(name='noise_factor_long', format='D', array=rms_lin_long)        

    fits_table = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7])    
    
    fits_table.writeto(rfi_table, overwrite = True)
  
    logger.info("\t ... RFI table saved ...\n")


def plot_rfi_im(cfg_par,time_step=-1):      
    '''
    
    Plots the .fits image output of rfi_im jpg format.
    
    '''        
    
    #check if image exists
    #plot image
    logger.info("\t ... Plotting RFI in 2D ... \n")

    if time_step != -1:
        time_tmp = int(float(cfg_par['rfi']['time_chunks']['time_step'])*time_step)
        if time_tmp == 0:
            time_name = '00'+str(time_tmp)+'m'
        elif time_tmp <100:
            time_name = '0'+str(time_tmp)+'m'
        else:
            time_name= str(time_tmp)+'m'    
    else:
        time_name = 'full' 

    rfi_freq_base = cfg_par['general']['rfidir']+'freq_base_'+time_name+'_im.fits'   
    rfi_freq_base_plot = cfg_par['general']['plotdir']+'freq_base_'+time_name+'_pl.png'
    


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
        time_tmp = int(float(cfg_par['rfi']['time_chunks']['time_step'])*time_step)
        if time_tmp == 0:
            time_name = '00'+str(time_tmp)+'m'
        elif time_tmp <100:
            time_name = '0'+str(time_tmp)+'m'
        else:
            time_name= str(time_tmp)+'m'    
    else:
        time_name = 'full' 

    rfi_freq_base = cfg_par['general']['rfidir']+'freq_base_'+time_name+'_im.fits'   
    rfi_freq_base_plot = cfg_par['general']['plotdir']+'freq_base_'+time_name+'_pl.png'
    

    # open file
    t = pyfits.open(rfi_freq_base)
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
    cbar.set_label(r'$\% > 5 \times$ r.m.s.',size=20)



    times, start, end = rfi.time_chunk(cfg_par)
    start.format='iso' 
    start.subformat='date_hm'   
    end.format='iso'
    end.subformat='date_hm'
    title_plot = '{0:%d}{0:%b}{0:%y}: {0:%H}:{0:%M} - {1:%H}:{1:%M}'.format(start.datetime,end.datetime)
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

    if time_step != -1:
        time_tmp = int(float(cfg_par['rfi']['time_chunks']['time_step'])*time_step)
        if time_tmp == 0:
            time_name = '00'+str(time_tmp)+'m'
        elif time_tmp <100:
            time_name = '0'+str(time_tmp)+'m'
        else:
            time_name= str(time_tmp)+'m'    
    else:
        time_name = 'full' 

    table_name = str(table_tmp[0])+'_'+time_name+'.fits'
    rfi_table = cfg_par['general']['tabledir']+table_name
    

    t = pyfits.open(rfi_table)
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
        rms = float(cfg_par['rfi']['theo_rms'])
        noise_all = noise_factor*rms
        noise_short = noise_factor_short*rms
        noise_long = noise_factor_long*rms

    if cfg_par['plots']['plot_noise'] == 'noise_factor':
        logger.info("\t ... Plotting factor of noise increas per frequency channel ...")
        noise_all = noise_factor
        noise_short = noise_factor_short
        noise_long = noise_factor_long    
        out_plot = cfg_par['general']['plotdir']+'freq_noise_'+time_name
    if cfg_par['plots']['plot_noise'] == 'flag':
        logger.info("\t ... Plotting percentage of flagged RFI per frequency channel ...")
        noise_all = flags
        noise_long = flags_long
        noise_short = flags_short
        out_plot = cfg_par['general']['plotdir']+'freq_flags_'+time_name


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
    
    if cfg_par['plots']['plot_noise'] != 'flag':
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
        ax1.set_ylabel(r'Factor of noise increase')
        ax1.get_yaxis().set_major_formatter(ticker.ScalarFormatter())

    if cfg_par['plots']['plot_noise'] == 'noise':
        ax1.set_yticks([1,2,3,5,10,50]) 
        ax1.set_ylabel(r'Predicted noise [mJy beam$^{-1}$]')     
        out_plot = out_plot+'_noise'+self.aperfi_plot_format    
        ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

    if cfg_par['plots']['plot_noise'] == 'flag':
        ax1.set_ylabel(r'$\% >$ '+str(cfg_par['rfi']['rms_clip'])+'*rms') 
    
    legend = plt.legend()
    legend.get_frame().set_edgecolor('black')

    times, start, end = rfi.time_chunk(cfg_par)
    start.format='iso' 
    start.subformat='date_hm'   
    end.format='iso'
    end.subformat='date_hm'
    title_plot = '{0:%d}{0:%b}{0:%y}: {0:%H}:{0:%M} - {1:%H}:{1:%M}'.format(start.datetime,end.datetime)
    ax1.set_title(title_plot)

    # Save figure to file
    rfi_freq_plot = out_plot+'_pl.png'
    plt.savefig(rfi_freq_plot,format='png',overwrite = True)      

    logger.info("\t ... RFI in 1D plotted ...\n\n")
   