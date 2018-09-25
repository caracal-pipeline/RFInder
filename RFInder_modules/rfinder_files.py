import os,string,sys
import numpy as np
from astropy.io import fits as fits
from astropy import units as u



import rfinder_stats as rfi_stats 
rfiST = rfi_stats.rfi_stats()

import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def write_freq_base(cfg_par,rms,time_step=-1) :
    '''
    
    Writes an image.fits of visibilities ordered by frequency, baseline.
    Baselines are ordered by their length.
    
    '''
    logger.info('\t ... Wrtiting RFI by frequency and baseline ...')


    #reverse frequencies if going from high-to-low         
    if time_step != -1:
        outputdir = cfg_par['general']['rfitimedir']
        time_tmp = int(float(cfg_par['rfi']['chunks']['time_step'])*time_step)
        if time_tmp == 0:
            time_name = '00'+str(time_tmp)+'m'
        elif time_tmp <100:
            time_name = '0'+str(time_tmp)+'m'
        else:
            time_name= str(time_tmp)+'m'    
    else:
        outputdir = str(cfg_par['general']['rfidir'])
        time_name = 'full' 

    if cfg_par['rfi']['RFInder_mode'] == 'use_flags':
        rfi_freq_base =outputdir+'flags_base_'+time_name+'.fits'
    if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
        rfi_freq_base =outputdir+'rfi_base_'+time_name+'.fits'

    #set fits file
    hdu = fits.PrimaryHDU(rms)
    hdulist = fits.HDUList([hdu])
    header = hdulist[0].header
    cfg_par['rfi']['chan_widths']
    #write header keywords               
    header['CRPIX1'] = 1
    header['CDELT1'] = np.abs(float(cfg_par['rfi']['chan_widths'])/1e6)
    header['CRVAL1'] = float(cfg_par['rfi']['lowfreq'])/1e6
    header['CTYPE1'] = ('Frequency')
    header['CUNIT1'] = ('MHz')
    header['CRPIX2'] = 1
    header['CRVAL2'] = 1
    header['CDELT2'] = 1
    header['CTYPE2'] =  ('Baseline')
    header['POLAR'] =  ('XX')
    header['BUNIT'] = ('% > '+str(5)+'*rms')
    header['BTYPE'] = ('intensity')
    #write file
    fits.writeto(rfi_freq_base,rms,header,overwrite=True)
    hdulist.close()

    logger.info('\t ... RFI on fits file ...\n')


    return 0


def rfi_frequency(cfg_par,time_step=-1):
    '''
    Determines the rfi per frequency channel. Saves results in table rfi_table.fits
    For each channel the flag % and factor of noise increase are stored for all, long and short baselines
    Long and short baselines are separated in half, depending on the number of baselines
    '''
    table_tmp = string.split(cfg_par['general']['msname'][0],'.MS')
    if len(table_tmp) == 1:
        table_tmp = string.split(cfg_par['general']['msname'][0],'.ms')

    if time_step != -1:
        tabledir = cfg_par['general']['timetabledir'] 
        outputdir = cfg_par['general']['rfitimedir'] 
        time_tmp = int(float(cfg_par['rfi']['chunks']['time_step'])*time_step)
        if time_tmp == 0:
            time_name = '00'+str(time_tmp)+'m'
        elif time_tmp <100:
            time_name = '0'+str(time_tmp)+'m'
        else:
            time_name= str(time_tmp)+'m'
    else:
        outputdir = cfg_par['general']['rfidir']
        tabledir = cfg_par['general']['tabledir']        
        time_name = 'full'
    
    if cfg_par['rfi']['RFInder_mode']== 'use_flags':
        rfi_freq_base =outputdir+'flags_base_'+time_name+'.fits'
        table_name = str(table_tmp[0])+'_flags_'+time_name+'.fits'
    if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
        rfi_freq_base =outputdir+'rfi_base_'+time_name+'.fits'
        table_name = str(table_tmp[0])+'_rfi_'+time_name+'.fits'
    
    rfi_table = tabledir+table_name

    #open file
    if os.path.exists(rfi_freq_base) == False:
        logger.error('### Image of RFI sorted by frequency over baseline lenght does not exist ###')    
    else:    
        
        # read data and header
        hdulist = fits.open(rfi_freq_base)  # read input                
        datacube = hdulist[0].data    
        prihdr = hdulist[0].header

        #set array of frequencies
        freqs = (np.linspace(1, datacube.shape[1], datacube.shape[1])\
                     - prihdr['CRPIX1'])*prihdr['CDELT1'] + prihdr['CRVAL1']
        
        # set y-array
        rms_lin = np.zeros([datacube.shape[1]])    
        flag_lin = np.zeros([datacube.shape[1]])   

        more_long = float(cfg_par['rfi']['number_baseline'])/float(cfg_par['rfi']['num_long'])
        more_short = float(cfg_par['rfi']['number_baseline'])/float(cfg_par['rfi']['num_short'])

        rms_lin_long = np.zeros([datacube.shape[1]]) + np.sqrt(more_long)          
        rms_lin_short = np.zeros([datacube.shape[1]]) + np.sqrt(more_short)

        flag_lin_long = np.zeros([datacube.shape[1]]) + (100.-100./more_long)    
        flag_lin_short = np.zeros([datacube.shape[1]]) + (100.-100./more_short)

        natural_rms = cfg_par['rfi']['theo_rms'] 

        elevation = np.zeros(freqs.shape)+cfg_par['rfi']['altaz'].alt*u.deg
        azimuth = np.zeros(freqs.shape)+cfg_par['rfi']['altaz'].az*u.deg


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
        if cfg_par['rfi']['chunks']['spw_enable'] == True:

            
            if cfg_par['rfi']['RFInder_mode']== 'use_flags':
                table_name_bin = str(table_tmp[0])+'_flags_'+time_name+'_spwbin.fits'
            if cfg_par['rfi']['RFInder_mode']== 'rms_clip':
                table_name_bin = str(table_tmp[0])+'_rfi_'+time_name+'_spwbin.fits'


            rfi_table_bin = tabledir+table_name_bin

            step_bin = cfg_par['rfi']['chunks']['spw_width']

            freqs_bin=np.arange(freqs[0],freqs[-1]+step_bin,step_bin)

            flag_lin_bin=np.zeros(freqs_bin.shape)
            rms_lin_bin=np.zeros(freqs_bin.shape)
            flag_lin_bin_short=np.zeros(freqs_bin.shape)
            rms_lin_bin_short=np.zeros(freqs_bin.shape)
            flag_lin_bin_long=np.zeros(freqs_bin.shape)
            rms_lin_bin_long=np.zeros(freqs_bin.shape)
            # is this correct ???
            natural_rms_bin=np.zeros(freqs_bin.shape)

            elevation_bin = np.zeros(freqs_bin.shape)+cfg_par['rfi']['altaz'].alt*u.deg
            azimuth_bin = np.zeros(freqs_bin.shape)+cfg_par['rfi']['altaz'].az*u.deg

            for i in xrange(0, len(freqs_bin)-1):
                #look for the right velocity bin
                index = (freqs_bin[i] <= freqs) & (freqs < freqs_bin[i+1])

                flag_lin_bin[i] = np.nanmean(flag_lin[index])
                rms_lin_bin[i] = np.nanmean(rms_lin[index])
                flag_lin_bin_short[i] = np.nanmean(flag_lin_short[index])
                rms_lin_bin_short[i] = np.nanmean(rms_lin_short[index])
                flag_lin_bin_long[i] = np.nanmean(flag_lin_long[index])
                rms_lin_bin_long[i] = np.nanmean(rms_lin_long[index])
                natural_rms_bin[i] = np.nanmean(natural_rms[index])

            # save fits table        
            c1 = fits.Column(name='frequency', format='D', unit='MHz', array=freqs_bin)
            c2 = fits.Column(name='percentage_flags', format='D', unit='-', array=flag_lin_bin)
            c3 = fits.Column(name='noise_factor', format='D', unit = '-', array=rms_lin_bin)
            c4 = fits.Column(name='percentage_flags_short', format='D', unit='-', array=flag_lin_bin_short)
            c5 = fits.Column(name='noise_factor_short', format='D', unit = '-', array=rms_lin_bin_short)
            c6 = fits.Column(name='percentage_flags_long', format='D', unit='-', array=flag_lin_bin_long)
            c7 = fits.Column(name='noise_factor_long', format='D', array=rms_lin_bin_long)        
            c8 = fits.Column(name='natural_rms', format='D', array=natural_rms_bin)        
            c9 = fits.Column(name='altitude', format='D', unit='deg', array=elevation_bin)        
            c10 = fits.Column(name='azimuth', format='D', unit='deg', array=azimuth_bin)        

            fits_table = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8, c9, c10])    

            fits_table.writeto(rfi_table_bin, overwrite = True)


        # save fits table        
        c1 = fits.Column(name='frequency', format='D', unit='MHz', array=freqs)
        c2 = fits.Column(name='percentage_flags', format='D', unit='-', array=flag_lin)
        c3 = fits.Column(name='noise_factor', format='D', unit = '-', array=rms_lin)
        c4 = fits.Column(name='percentage_flags_short', format='D', unit='-', array=flag_lin_short)
        c5 = fits.Column(name='noise_factor_short', format='D', unit = '-', array=rms_lin_short)
        c6 = fits.Column(name='percentage_flags_long', format='D', unit='-', array=flag_lin_long)
        c7 = fits.Column(name='noise_factor_long', format='D', array=rms_lin_long)        
        c8 = fits.Column(name='natural_rms', format='D', array=natural_rms)        
        c9 = fits.Column(name='altitude', format='D', unit='deg', array=elevation)        
        c10 = fits.Column(name='azimuth', format='D', unit='deg', array=azimuth)        
        
        fits_table = fits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7, c8])    
        
        fits_table.writeto(rfi_table, overwrite = True)
      
        logger.info("\t ... RFI table saved ...\n")


