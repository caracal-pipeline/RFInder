import logger
from matplotlib import gridspec
from matplotlib import pyplot as plt
import aplpy


def rfi_frequency(cfg_par):
    '''
    Determines the rfi per frequency channel. Saves results in table rfi_table.fits
    For each channel the flag % and factor of noise increase are stored for all, long and short baselines
    Long and short baselines are separated in half, depending on the number of baselines
    '''

    #open file
    if os.path.exists(self.rfi_freq_base) == False:
        self.logger.error('### Image of RFI sorted by frequency over baseline lenght does not exist ###')    
        self.logger.error('### Run aperfi.rfi_im() first ###')  
    else:    
        
        # read data and header
        hdulist = pyfits.open(self.rfi_freq_base)  # read input                
        datacube = hdulist[0].data    
        prihdr = hdulist[0].header

        #set array of frequencies
        self.freqs = (np.linspace(1, datacube.shape[1], datacube.shape[1])\
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

            shortbase=datacube[:int(datacube.shape[0]/2),i]
            longbase = datacube[int(datacube.shape[0]/2):,i]               
            
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
    
        # save fits table        
        c1 = pyfits.Column(name='frequency', format='D', unit='MHz', array=self.freqs)
        c2 = pyfits.Column(name='flag', format='D', unit='-', array=flag_lin)
        c3 = pyfits.Column(name='noise_factor', format='D', unit = '-', array=rms_lin)
        c4 = pyfits.Column(name='flag_short', format='D', unit='-', array=flag_lin_short)
        c5 = pyfits.Column(name='noise_factor_short', format='D', unit = '-', array=rms_lin_short)
        c6 = pyfits.Column(name='flag_long', format='D', unit='-', array=flag_lin_long)
        c7 = pyfits.Column(name='noise_factor_long', format='D', array=rms_lin_long)        

        fits_table = pyfits.BinTableHDU.from_columns([c1, c2, c3, c4, c5, c6, c7])    
        
        fits_table.writeto(self.rfi_table, overwrite = True)
        

def plot_rfi_im(self):      
    '''
    
    Plots the .fits image output of rfi_im jpg format.
    
    '''        
    
    #check if image exists
    #plot image
    fig = aplpy.FITSFigure(self.rfi_freq_base,figsize=(12,8))

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
    plt.savefig(self.rfi_freq_base_plot,format='png' ,overwrite=True)


def plot_noise_frequency(self):
    '''
    Plots the noise or % of rfi per frequency channel for all, long and short baselines.
    In default.cfga
        aperfi_noise = 'rfi' (or flag)
        aperfi_plot_long_short = False
    '''

    #open file
    if os.path.exists(self.rfi_table) == False:
        self.logger.error('### Table of RFI and flags of visibilities does not exist ###')    
        self.logger.error('### Run aperfi.rfi_frequency() first ###')  
    else:  


        t = pyfits.open(self.rfi_table)
        data_vec = t[1].data
        cols = t[1].columns
        
        freqs = np.array(data_vec['frequency'],dtype=float)
        flags = np.array(data_vec['flag'],dtype=float)
        noise_factor = np.array(data_vec['noise_factor'],dtype=float)
        noise_factor_long = np.array(data_vec['noise_factor_long'],dtype=float)
        flags_long = np.array(data_vec['flag_long'],dtype=float)
        noise_factor_short = np.array(data_vec['noise_factor_short'],dtype=float)
        flags_short = np.array(data_vec['flag_short'],dtype=float)

       
        #if self.aperfi_noise == 'noise':
        #    self.predicted_noise_channel()
        #    noise_all = noise_factor*self.noise_freq
        #    noise_short = noise_factor_short*self.noise_freq
        #    noise_long = noise_factor_long*self.noise_freq
        if self.aperfi_plot_noise == 'rfi':
            noise_all = noise_factor
            noise_short = noise_factor_short
            noise_long = noise_factor_long          
        if self.aperfi_plot_noise == 'flag':
            noise_all = flags
            noise_long = flags_long
            noise_short = flags_short


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
        
        if self.aperfi_noise != 'flag':
            ax1.set_yscale('log', basey=10)
    
        #define title output                     
        
        #plot
        label_all = 'All baselines' 
        label_long = 'Long baselines' 
        label_short = 'Short baselines' 

        if self.aperfi_plot_long_short == True:
            ax1.step(freqs,noise_short, where= 'pre', color='red', linestyle='-',label=label_short)
            ax1.step(freqs,noise_long, where= 'pre', color='blue', linestyle='-',label=label_long)
            out_plot = out_plot+'_sl_'

        ax1.step(freqs,noise_all, where= 'pre', color='black', linestyle='-',label=label_all)

        #titleplot = self.target+': '+self.aperfi_startime+' - '+self.aperfi_endtime
        #ax1.set_title(titleplot)
        
        # set axis, legend ticks

        ax1.set_xlim([np.min(freqs)-5,np.max(freqs)+5])
        xticks_num = np.linspace(int(self.channelFreqs[0]),int(self.channelFreqs[-1]),10,dtype=int)
        ax1.set_xticks(xticks_num)

        if self.aperfi_plot_noise == 'rfi':
            ax1.set_yticks([1,round(np.sqrt(2),2),2,3,5,10,50]) 
            ax1.set_ylabel(r'Factor of noise increase')
            ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        # if self.aperfi_noise == 'noise':
        #     ax1.set_yticks([1,2,3,5,10,50]) 
        #     ax1.set_ylabel(r'Predicted noise [mJy beam$^{-1}$]')     
        #     out_plot = out_plot+'_noise'+self.aperfi_plot_format    
        #     ax1.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

        if self.aperfi_plot_noise == 'flag':
            ax1.set_ylabel(r'$\% >$ '+str(self.aperfi_rmsclip)+'*rms') 
        
        legend = plt.legend()
        legend.get_frame().set_edgecolor('black')

        # Save figure to file
        plt.savefig(self.rfi_freq_plot,format='png',overwrite = True)