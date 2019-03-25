import os, string, sys, glob
import numpy as np
from astropy.io import fits as fits
from astropy import units as u
from jinja2 import FileSystemLoader, Environment
import shutil
import rfinder_stats as rfi_stats 
rfiST = rfi_stats.rfi_stats()

import logging
logger = logging.getLogger('log-rfinder.log')
#logger.setLevel(logging.INFO)

#fh = logging.FileHandler('log-rfinder.log')
#fh.setLevel(logging.INFO)

#ch = logging.StreamHandler()
#ch.setLevel(logging.WARNING)

#formatter = logging.Formatter('%(levelname)s - %(filename)s - %(message)s')
#fh.setFormatter(formatter)
#ch.setFormatter(formatter)

#logger.addHandler(ch)
#logger.addHandler(fh)


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

        elevation = np.zeros(freqs.shape)+cfg_par['rfi']['altaz'].alt
        azimuth = np.zeros(freqs.shape)+cfg_par['rfi']['altaz'].az


        for i in xrange(0,datacube.shape[1]):
            
            flag_lin_tmp = np.divide(np.sum(datacube[:,i]),datacube.shape[0])
            flag_lin[i] = flag_lin_tmp


            baseline_cutoff = float(cfg_par['rfi']['baseline_cut'])
            lenghts = np.array([cfg_par['rfi']['baseline_lenghts']])+0.
            idx = (np.abs(lenghts - baseline_cutoff)).argmin()

            shortbase=datacube[:idx,i]
            longbase = datacube[idx:,i]               
            
            rms_lin_tmp = 1.-np.divide(np.divide(np.nansum(datacube[:,i]),datacube.shape[0]),100.)
            rms_lin[i] = np.divide(1.,np.sqrt(rms_lin_tmp))

            flag_lin_tmp = np.divide(np.nansum(shortbase),len(shortbase))
            flag_lin_short[i] = flag_lin_tmp
            rms_lin_tmp_short = 1.-np.divide(np.divide(np.nansum(shortbase),len(shortbase)),100.)
            rms_lin_short[i] *= np.divide(1.,np.sqrt(rms_lin_tmp_short))

            flag_lin_tmp = np.divide(np.nansum(longbase),len(longbase))
            flag_lin_long[i] = flag_lin_tmp
            rms_lin_tmp_long = 1.-np.divide(np.divide(np.nansum(longbase),len(longbase)),100.)
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

        return 0

# Content to be published
def write_html_fullreport(cfg_par):

    # Configure Jinja and ready the template
    env = Environment(
        loader=FileSystemLoader(cfg_par['general']['template_folder'])
    )

    #base_template = env.get_template('report.html')
    # Content to be published
    title = 'RFI report: {0:s}'.format(cfg_par['general']['msname'])

    #imagename1 = '/Users/maccagni/Projects/RFI/rfinder_test/rfi/plots/altaz/AltAZ_rfi1297-1317MHz.png'
    #data_uri1 = open(imagename1, 'rb').read().encode('base64').replace('\n', '')
    
    if cfg_par['rfi']['RFInder_mode'] == 'rms_clip':
        imagename1 = cfg_par['general']['plotdir']+'rfi_base_full.png'
        imagename3 = cfg_par['general']['plotdir']+'noise_full_sl_rfi.png'
    elif cfg_par['rfi']['RFInder_mode'] == 'use_flags':
        imagename1 = cfg_par['general']['plotdir']+'flags_base_full.png'
        imagename3 = cfg_par['general']['plotdir']+'noise_full_sl_flags.png'
    
    data_uri1 = open(imagename1, 'rb').read().encode('base64').replace('\n', '')

    imagename2 = cfg_par['general']['plotdir']+'AltAZ_full.png'
    data_uri2 = open(imagename2, 'rb').read().encode('base64').replace('\n', '')

    data_uri3 = open(imagename3, 'rb').read().encode('base64').replace('\n', '')

    video_name1 = cfg_par['general']['moviedir']+'AltAz_movie.gif'
    if os.path.exists(video_name1):
        video_encoded1 = open(video_name1, "rb").read().encode("base64")

    video_name2 = cfg_par['general']['moviedir']+'Time_2Dplot_movie.gif'
    if os.path.exists(video_name2):
        video_encoded2 = open(video_name2, "rb").read().encode("base64")

    video_name3 = cfg_par['general']['moviedir']+'TimeChunks_1D_noise.gif'
    if os.path.exists(video_name3):
        video_encoded3 = open(video_name3, "rb").read().encode("base64")

    if cfg_par['plots']['movies']['movies_in_report'] == True:
        template = env.get_template('full_template.html')
        with open(cfg_par['general']['rfidir']+'full_report.html', "w") as f:
            lenghts = np.array([cfg_par['rfi']['baseline_lenghts']])+0.
            f.write(template.render(
                title=title,
                fieldname=cfg_par['general']['fieldname'],
                field=cfg_par['general']['field'],
                totchans = int(cfg_par['rfi']['total_channels']),
                chan_widths=round(cfg_par['rfi']['chan_widths']/1e3,4),
                lowfreq=round(cfg_par['rfi']['lowfreq']/1e6,3),
                highfreq=round(cfg_par['rfi']['highfreq']/1e6,3),
                startdate = ('{0:%y}{0:%b}{0:%d} {0:%X}'.format(cfg_par['rfi']['startdate'].datetime)),
                enddate =   ('{0:%y}{0:%b}{0:%d} {0:%X}'.format(cfg_par['rfi']['enddate'].datetime)),
                nant = cfg_par['rfi']['nant'],
                ant_names = cfg_par['rfi']['ant_names'],
                maxbase = str(np.round(lenghts[0][-1],0)),
                minbase = str(np.round(lenghts[0][0],0)),
                totbase = cfg_par['rfi']['number_baseline'],
                exptime = np.round(cfg_par['rfi']['exptime'],2),
                polnum = cfg_par['rfi']['polnum'],
                noise = np.round(cfg_par['rfi']['theo_rms'][0]*1e3,5),
                img_tag1 = '<img class="a" src="data:image/png;base64,{0}">'.format(data_uri1),
                video_tag1 = '<img class="b" src="data:video/gif;base64,{0}">'.format(video_encoded1),
                img_tag3 = '<img class="c" src="data:image/png;base64,{0}">'.format(data_uri3),
                #video_tag1 = '<img class="d" src="data:video/gif;base64,{0}">'.format(video_encoded1),        
                video_tag2 = '<img class="e" src="data:video/gif;base64,{0}">'.format(video_encoded2),        
                video_tag3 = '<img class="f" src="data:video/gif;base64,{0}">'.format(video_encoded3)                    
            ))

    elif cfg_par['plots']['movies']['movies_in_report'] == False:
        template = env.get_template('fullshort_template.html')
        with open(cfg_par['general']['rfidir']+'full_report.html', "w") as f:
            lenghts = np.array([cfg_par['rfi']['baseline_lenghts']])+0.
            f.write(template.render(
                title=title,
                fieldname=cfg_par['general']['fieldname'],
                field=cfg_par['general']['field'],
                totchans = int(cfg_par['rfi']['total_channels']),
                chan_widths=round(cfg_par['rfi']['chan_widths']/1e3,4),
                lowfreq=round(cfg_par['rfi']['lowfreq']/1e6,3),
                highfreq=round(cfg_par['rfi']['highfreq']/1e6,3),
                startdate = ('{0:%y}{0:%b}{0:%d} {0:%X}'.format(cfg_par['rfi']['startdate'].datetime)),
                enddate =   ('{0:%y}{0:%b}{0:%d} {0:%X}'.format(cfg_par['rfi']['enddate'].datetime)),
                nant = cfg_par['rfi']['nant'],
                ant_names = cfg_par['rfi']['ant_names'],
                maxbase = str(np.round(lenghts[0][-1],0)),
                minbase = str(np.round(lenghts[0][0],0)),
                totbase = cfg_par['rfi']['number_baseline'],
                exptime = np.round(cfg_par['rfi']['exptime'],2),
                polnum = cfg_par['rfi']['polnum'],
                noise = np.round(cfg_par['rfi']['theo_rms'][0]*1e3,5),
                img_tag1 = '<img class="a" src="data:image/png;base64,{0}">'.format(data_uri1),
                video_tag1 = '<img class="b" src="data:video/gif;base64,{0}">'.format(data_uri2),
                img_tag3 = '<img class="c" src="data:image/png;base64,{0}">'.format(data_uri3),                 
            ))

        logger.info('---- Html report done ----\n')

        return 0

def write_html_timereport(cfg_par):

    # Configure Jinja and ready the template
    env = Environment(
        loader=FileSystemLoader(cfg_par['general']['template_folder'])
    )
    template = env.get_template('time_template.html')

    #base_template = env.get_template('report.html')
    # Content to be published
    title = 'RFI time scans report: {0:s}'.format(cfg_par['general']['msname'])

    #imagename1 = '/Users/maccagni/Projects/RFI/rfinder_test/rfi/plots/altaz/AltAZ_rfi1297-1317MHz.png'
    #data_uri1 = open(imagename1, 'rb').read().encode('base64').replace('\n', '')
    video_name1 = cfg_par['general']['moviedir']+'AltAz_movie.gif'
    if os.path.exists(video_name1):
        video_encoded1 = open(video_name1, "rb").read().encode("base64")
    else:
        video_encoded1 = None
    
    video_name2 = cfg_par['general']['moviedir']+'Time_2Dplot_movie.gif'
    if os.path.exists(video_name2):
        video_encoded2 = open(video_name2, "rb").read().encode("base64")
    else:
        video_encoded1 = None

    video_name3 = cfg_par['general']['moviedir']+'TimeChunks_1D_noise.gif'
    if os.path.exists(video_name3):
        video_encoded3 = open(video_name3, "rb").read().encode("base64")
    else:
        video_encoded3 = None

    if cfg_par['plots']['movies']['movies_in_report'] == True:

        with open(cfg_par['general']['rfidir']+'time_report.html', "w") as f:
            lenghts = np.array([cfg_par['rfi']['baseline_lenghts']])+0.
            f.write(template.render(
                title=title,
                fieldname=cfg_par['general']['fieldname'],
                field=cfg_par['general']['field'],
                totchans = int(cfg_par['rfi']['total_channels']),
                chan_widths=round(cfg_par['rfi']['chan_widths']/1e3,4),
                lowfreq=round(cfg_par['rfi']['lowfreq']/1e6,3),
                highfreq=round(cfg_par['rfi']['highfreq']/1e6,3),
                startdate = ('{0:%y}{0:%b}{0:%d} {0:%X}'.format(cfg_par['rfi']['startdate'].datetime)),
                enddate =   ('{0:%y}{0:%b}{0:%d} {0:%X}'.format(cfg_par['rfi']['enddate'].datetime)),
                nant = cfg_par['rfi']['nant'],
                ant_names = cfg_par['rfi']['ant_names'],
                maxbase = str(np.round(lenghts[0][-1],0)),
                minbase = str(np.round(lenghts[0][0],0)),
                totbase = cfg_par['rfi']['number_baseline'],
                exptime = np.round(cfg_par['rfi']['exptime']*60.,2),
                polnum = cfg_par['rfi']['polnum'],
                noise = np.round(cfg_par['rfi']['theo_rms'][0]*1e3,5),
                video_tag1 = '<img class="b" src="data:video/gif;base64,{0}">'.format(video_encoded1),        
                video_tag2 = '<img class="a" src="data:video/gif;base64,{0}">'.format(video_encoded2),        
                video_tag3 = '<img class="c" src="data:video/gif;base64,{0}">'.format(video_encoded3)                    
            ))

    elif cfg_par['plots']['movies']['movies_in_report'] == False:

        logger.info('\t ERROR:  movies in report must be set to TRUE\n')

        return 0

    logger.info('---- Html time report done ----\n')

    return 0

def find_altaz_plots(cfg_par):
    
    tmp_arr=[]

    if cfg_par['rfi']['RFInder_mode']=='use_flags':
        filenames = glob.glob(cfg_par['general']['altazplotdir']+'/AltAZ_flags*')

        for i in xrange(0,len(filenames)):
            tmp = filenames[i].split('_flags')
            tmp_arr.append(tmp[1].split('MHz')[0])
        tmp_arr.sort()
            
        filenames = [tmp[0]+'_flags' + s for s in tmp_arr] 
        filenames = [s + 'MHz.png' for s in filenames] 
    
    elif cfg_par['rfi']['RFInder_mode']=='rms_clip':
        filenames = glob.glob(cfg_par['general']['altazplotdir']+'/AltAZ_rfi*')

        for i in xrange(0,len(filenames)):
            tmp = filenames[i].split('_rfi')
            tmp_arr.append(tmp[1].split('MHz')[0])
        tmp_arr.sort()
            
        filenames = [tmp[0]+'_rfi' + s for s in tmp_arr] 
        filenames = [s + 'MHz.png' for s in filenames] 

    return filenames

def find_2d_plots(cfg_par):

    if cfg_par['rfi']['RFInder_mode']=='use_flags':
        filenames = glob.glob(cfg_par['general']['timeplotdir2D']+'/flags_base*')
    elif cfg_par['rfi']['RFInder_mode']=='rms_clip':
        filenames = glob.glob(cfg_par['general']['timeplotdir2D']+'/rfi_base*')

    tmp_arr=[]
    for i in xrange(0,len(filenames)):
            tmp = filenames[i].split('base_')[1]
            tmp_arr.append(tmp.split('m.png')[0])
    tmp_arr.sort()

    tmp = filenames[0].split('base_')
    filenames = [tmp[0]+'base_' + s for s in tmp_arr] 
    filenames = [s + 'm.png' for s in filenames] 

    return filenames

def find_1d_plots(cfg_par,name_root):

    #select files
    filenames = glob.glob(cfg_par['general']['timeplotdir1D']+'/'+name_root+'_*')
    #print filenames 

    tmp_arr=[]

    for i in xrange(0,len(filenames)):
        if cfg_par['rfi']['RFInder_mode']=='use_flags':
            if len(filenames[i].split('sl_rfi.'))>1:
                continue
        elif cfg_par['rfi']['RFInder_mode']=='rms_clip':
            if  len(filenames[i].split('sl_flags.'))>1:
                continue
        tmp = filenames[i].split(name_root+'_')
        tmp_arr.append(tmp[1].split('m_sl_rfi.png')[0])

    tmp_arr.sort()
    tmp = filenames[0].split(name_root+'_')

    filenames = [tmp[0]+name_root+'_' + s for s in tmp_arr]
    #if cfg_par['rfi']['RFInder_mode']=='use_flags':
    #    filenames = [s + 'm_sl_flags.png' for s in filenames] 
    #elif cfg_par['rfi']['RFInder_mode']=='rms_clip':
    #    filenames = [s + 'm_sl_rfi.png' for s in filenames] 
    
    #print filenames

    return filenames

def set_dirs(cfg_par):
    '''
 
    Sets temporary subtirectories
    Creates directory rfi_stokes/
    Subdirectories can be kept using no_cleanup option

    '''

    key = 'general'

    if os.path.exists(cfg_par[key]['moviedir']) == False:
         os.makedirs(cfg_par[key]['moviedir'])

    if os.path.exists(cfg_par[key]['rfidir']) == False:
         os.makedirs(cfg_par[key]['rfidir'])           

    if os.path.exists(cfg_par[key]['tabledir']) == False:
         os.makedirs(cfg_par[key]['tabledir'])

    if os.path.exists(cfg_par[key]['plotdir']) == False:
         os.makedirs(cfg_par[key]['plotdir'])

    if cfg_par['rfi']['chunks']['time_enable'] == True:

        if os.path.exists(cfg_par[key]['rfitimedir']) == False:
             os.makedirs(cfg_par[key]['rfitimedir'])

        if os.path.exists(cfg_par[key]['timetabledir']) == False:
             os.makedirs(cfg_par[key]['timetabledir'])

        if os.path.exists(cfg_par[key]['timechunksdir']) == False:
             os.makedirs(cfg_par[key]['timechunksdir'])

        if os.path.exists(cfg_par[key]['timeplotdir1D']) == False:
             os.makedirs(cfg_par[key]['timeplotdir1D'])

        if os.path.exists(cfg_par[key]['timeplotdir2D']) == False:
             os.makedirs(cfg_par[key]['timeplotdir2D'])

        if os.path.exists(cfg_par[key]['altazplotdir']) == False:
             os.makedirs(cfg_par[key]['altazplotdir'])

        return(0)

def cleanup(cfg_par):

    key = 'general'


    if cfg_par['rfi']['chunks']['time_enable'] == True:
        if os.path.exists(cfg_par[key]['timeplotdir1D']):
            shutil.rmtree(cfg_par[key]['timeplotdir1D'])

        if os.path.exists(cfg_par[key]['timeplotdir2D']):
            shutil.rmtree(cfg_par[key]['timeplotdir2D']) 

        if os.path.exists(cfg_par[key]['altazplotdir']):
            shutil.rmtree(cfg_par[key]['altazplotdir']) 

        if os.path.exists(cfg_par[key]['timechunksdir']):
            shutil.rmtree(cfg_par[key]['timechunksdir']) 

        if os.path.exists(cfg_par[key]['tabledir']):
            shutil.rmtree(cfg_par[key]['tabledir'])

        if os.path.exists(cfg_par[key]['rfitimedir']):
            shutil.rmtree(cfg_par[key]['rfitimedir'])

    else:
        if os.path.exists(cfg_par[key]['plotdir']):
            shutil.rmtree(cfg_par[key]['plotdir'])
        if os.path.exists(cfg_par[key]['tabledir']):
            shutil.rmtree(cfg_par[key]['tabledir'])
            
    fitsname = glob.glob(cfg_par[key]['rfidir']+'*.fits')

    if len(fitsname)>0:
        for i in xrange(0,len(fitsname)):
            os.remove(fitsname[i])

