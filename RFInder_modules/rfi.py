# Import modules
import sys, string, os
import numpy as np
import yaml
import json
import glob
#import pyrap.tables as tables
import casacore.tables as tables
import logging


from astropy.time import Time
import numpy as np
from astropy import units as u

from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn
import rfinder

import rfinder_stats as rfi_stats

rfiST = rfi_stats.rfi_stats()



__author__ = "Filippo Maccagni, Tom Oosterloo"
__copyright__ = "RadioLife"
__version__ = "1.0.0"
__email__ = "filippo.maccagni@gmail.com"
__status__ = "Development"

class rfi:

    def __init__(self):

        #set logger
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

        #self.logger.info("\t ... Initializing rfi.py ... \n\n")

    def time_chunk(self,cfg_par):

        self.logger.info("\t ...  Observing time Info ... \n")
        self.msfile = cfg_par['general']['msfullpath']
        
        t=tables.table(self.msfile)
        self.time = t.getcol('TIME')
        t.close()

        starttime= self.time[0]
        endtime=self.time[-1]
        time_chunk = float(cfg_par['rfi']['time_chunks']['time_step'])*60.
        times=np.arange(starttime,endtime+time_chunk*1,time_chunk)
        
        startdate=Time(starttime/3600./24.,format='mjd',scale='utc')
        startdate.format='iso' 
        startdate.subformat='date_hm'       

        enddate=Time(endtime/3600./24.,format='mjd',scale='utc')
        enddate.format='iso'        
        enddate.subformat='date_hm'       


        self.logger.info('\t Start date: {0:%y}{0:%b}{0:%d}:{0:%X}'.format(startdate.datetime))
        self.logger.info('\t End date  : {0:%y}{0:%b}{0:%d}:{0:%X} \n\n'.format(enddate.datetime))

        return times,startdate,enddate

    def load_from_ms(self,cfg_par,times=-1):
        '''

        Loads important columns from MS file
        From MS: 
        Field_ID, Data, Flag, Antenna1, Antenna2,
        From MS/ANTENNA:
        Position,Name
        From MS/SPECTRAL_WINDOW
        Chan_width, Chan_freq

        '''

        self.logger.info("\t ... Loading MSfile ...\n")
        
        self.msfile = cfg_par['general']['msfullpath']
        self.aperfi_badant = cfg_par['rfi']['bad_antenna'] 
        self.selectFieldID = cfg_par['general']['field']
        antennas = tables.table(self.msfile +'/ANTENNA')
        self.ant_pos = np.array(antennas.getcol('POSITION'))
        self.ant_wsrtnames = np.array(antennas.getcol('NAME'))
        
        self.ant_names = np.arange(0,self.ant_wsrtnames.shape[0],1)
        self.nant = len(self.ant_names)

        #logging
        self.logger.info("\tTotal number of antennas:\t"+str(self.nant))
        self.logger.info("\tAntenna names:\t\t"+str(self.ant_names))

        antennas.close()

        spw=tables.table(self.msfile+'/SPECTRAL_WINDOW')
        self.channelWidths=spw.getcol('CHAN_WIDTH')
        self.channelFreqs=spw.getcol('CHAN_FREQ')
        cfg_par['rfi']['chan_widths'] = self.channelWidths[0][0]
        cfg_par['rfi']['lowfreq'] = float(self.channelFreqs[0][0])
        cfg_par['rfi']['highfreq'] = float(self.channelFreqs[-1][-1])

        spw.close()
       
        self.logger.info("\tBandwidth [kHz]:\t"+str(cfg_par['rfi']['chan_widths']/1e3))
        self.logger.info("\tStart     [GHz]:\t"+str(cfg_par['rfi']['lowfreq']/1e9))
        self.logger.info("\tEnd       [GHz]:\t"+str(cfg_par['rfi']['highfreq']/1e9)+'\n')

        t=tables.table(self.msfile)


        if times !=0:
            value_end = times[1]
            value_start = times[0]

            t2 = tables.taql('select from $t where TIME < $value_end and TIME>$value_start')
            self.fieldIDs=t2.getcol('FIELD_ID')
            self.ant1=t2.getcol('ANTENNA1')
            self.ant2=t2.getcol('ANTENNA2')

            selection=self.fieldIDs==self.selectFieldID
            selection*=self.ant1!=self.ant2

            self.vis = t2.getcol('DATA')[selection]
            self.flag = t2.getcol('FLAG')[selection]
            self.interval = t2.getcol('INTERVAL')[selection]
            t2.close()
        
        else:
            self.fieldIDs=t.getcol('FIELD_ID')
            self.ant1=t.getcol('ANTENNA1')
            self.ant2=t.getcol('ANTENNA2')

            selection=self.fieldIDs==self.selectFieldID
            selection*=self.ant1!=self.ant2

            self.vis = t.getcol('DATA')[selection]
            self.interval = t.getcol('INTERVAL')[selection]
            self.flag = t.getcol('FLAG')[selection]

        t.close()

        #determine number of baselines
        nrAnt=np.unique(np.concatenate((self.ant1,self.ant2))).shape[0]
        
        if not self.aperfi_badant:
            nrbadant =len(int(self.aperfi_badant))
        else:
            nrbadant = 0.

        
        nrBaseline=(nrAnt-nrbadant)*(nrAnt-nrbadant-1)/2        
        cfg_par['rfi']['number_baseline'] = nrBaseline

        rfiST.predict_noise(cfg_par,self.channelWidths,self.interval,self.flag)

        self.logger.info("\t ... info from MS file loaded  \n\n")

        return 0

    def baselines_from_ms(self,cfg_par):
        '''
        Reads which baselines were used in the observations
        Stores them sorted by lenght in the array baselines_sort
        Creates the Matrix to analize visibilites on each baseline separately
        '''

        self.logger.info('\t ... Sorting baselines ... \n')


        #sort baselines by length
        baselines = []
        for i in range(0,self.nant-1) :
            for j in range(i+1,self.nant) :
                # Exclude anticorrelations
                if self.aperfi_badant != None:
                    if self.ant_names[i]!=self.ant_names[j] and any(x != self.ant_names[i] for x in self.aperfi_badant) and any(x != self.ant_names[i] for x in self.aperfi_badant):
                        #distances are in units of meters
                        xdiff = (self.ant_pos[i,0]-self.ant_pos[j,0])*(self.ant_pos[i,0]-self.ant_pos[j,0])
                        ydiff = (self.ant_pos[i,1]-self.ant_pos[j,1])*(self.ant_pos[i,1]-self.ant_pos[j,1])
                        zdiff = (self.ant_pos[i,2]-self.ant_pos[j,2])*(self.ant_pos[i,2]-self.ant_pos[j,2])

                        dist = np.sqrt(xdiff+ydiff+zdiff)
                        baselines.append(([self.ant_names[i],self.ant_names[j]],dist))
                else: 
                    if self.ant_names[i]!=self.ant_names[j]:
                        #distances are in units of meters
                        xdiff = (self.ant_pos[i,0]-self.ant_pos[j,0])*(self.ant_pos[i,0]-self.ant_pos[j,0])
                        ydiff = (self.ant_pos[i,1]-self.ant_pos[j,1])*(self.ant_pos[i,1]-self.ant_pos[j,1])
                        zdiff = (self.ant_pos[i,2]-self.ant_pos[j,2])*(self.ant_pos[i,2]-self.ant_pos[j,2])

                        dist = np.sqrt(xdiff+ydiff+zdiff)
                        baselines.append(([self.ant_names[i],self.ant_names[j]],dist))


        self.baselines_sort = sorted(baselines, key=lambda baselines: baselines[1])  

        cfg_par['rfi']['baseline_lenghts'] = np.array(self.baselines_sort)[:,1]

        # Define matrix of indecese of baselines                                         
        self.blMatrix=np.zeros((self.nant,self.nant),dtype=int)
        for i in range(0,len(self.baselines_sort)) :

            ant1 = self.baselines_sort[i][0][0]
            ant2 = self.baselines_sort[i][0][1]
            self.blMatrix[ant1,ant2] = i
            self.blMatrix[ant2,ant1] = i


        self.logger.info('\t ... Matrix of baselines sorted by length loaded ...\n\n')


    def priors_flag(self,cfg_par):
        '''

        Flags YY,XY,YX polarizations
        Flags autocorrelations
        Flags bad antennas set by aperfi_badant = [ x, y, z ]
        Stores information on baselines sorted by lenght in baselines_sort 
        Stores visibilities in a cube (baselines[sorted by length], frequency, time_scans of each baseline)

        '''

        self.logger.info('\t ... Flagging a-prioris  ...\n')


        self.datacube = np.zeros([len(self.baselines_sort),self.vis.shape[1],self.vis.shape[0]/(len(self.baselines_sort))])
        baseline_counter = np.zeros((self.nant,self.nant),dtype=int)
        #flag unused polarizations
        pol = cfg_par['rfi']['polarization']
        if (pol == 'xx' or pol == 'XX'):
            self.flag[:,:,1] = True #YY
            #self.flag[:,:,2] = True #XY
            #self.flag[:,:,3] = True #YX
        elif (pol == 'yy' or pol == 'YY'):
            self.flag[:,:,0] = True #YY
            self.flag[:,:,2] = True #XY
            self.flag[:,:,3] = True #YX
        elif (pol == 'xy' or pol == 'XY'):
            self.flag[:,:,0] = True #YY
            self.flag[:,:,1] = True #XY
            self.flag[:,:,3] = True #YX      
        elif (pol == 'yx' or pol == 'YX'):
            self.flag[:,:,0] = True #YY
            self.flag[:,:,1] = True #XY
            self.flag[:,:,2] = True #YX  
        elif (pol == 'q' or pol == 'QQ'):
            self.flag[:,:,2] = True #XY
            self.flag[:,:,3] = True #YX


        #flag autocorrelations and bad antennas
        for i in xrange(0,self.vis.shape[0]):
            

            if self.aperfi_badant != None:
                if (any(x == self.ant1[i] for x in self.aperfi_badant) or any(x == self.ant2[i] for x in self.aperfi_badant)):
                    self.flag[i,:,0] = True
            if self.ant1[i] == self.ant2[i]:
                self.flag[i,:,0] = True            
            else:
                a1 = self.ant1[i]
                a2 = self.ant2[i]
                indice=self.blMatrix[a1,a2]
                # Count how many visibilities per baseline
                counter=baseline_counter[a1,a2]
                # Put amplitude of visibility
                # In the right place in the new array
                if (pol == 'xx' or pol == 'XX'):
                    self.datacube[indice,:,counter]=np.abs(self.vis[i,:,0])
                elif (pol == 'yy' or pol == 'YY'):
                    self.datacube[indice,:,counter]=np.abs(self.vis[i,:,1])
                elif (pol == 'xy' or pol == 'XY'):
                    self.datacube[indice,:,counter]=np.abs(self.vis[i,:,2])
                elif (pol == 'yx' or pol == 'YX'):
                    self.datacube[indice,:,counter]=np.abs(self.vis[i,:,3])
                #elif (pol == 'q' or pol == 'Q'):
                #    self.datacube[indice,:,counter]=np.abs(self.vis[i,:,0])-np.abs(self.vis[i,:,1])/2.

                # Update the number of visibility in that baseline
                baseline_counter[a1,a2]+=1

        #why is this here?
        #self.datacube = np.transpose(self.datacube, (1, 0, 2)) 

 

        self.logger.info('\t ... Prior flagging done ...\n\n')

        return self.datacube

    def find_rfi(self,datas,cfg_par,time_step=-1):
        '''

        For each baseline finds all signal above rms*aperfi_clip
        Creates a cube of visibilities sorted by baseline_lengh , frequency, time
        Stores them sorted by lenght in the array baselines_sort
        Creates the Matrix to analize visibilites on each baseline separately

        '''
        
        self.logger.info('\t ... Searching for RFI ...\n')
       
        self.datacube = datas

        rms = np.zeros([self.datacube.shape[0],self.datacube.shape[1]])
        self.mean_array = np.zeros([self.datacube.shape[0],self.datacube.shape[1]])
        self.flag_lim_array= np.zeros([self.datacube.shape[0]])

        self.aperfi_rfifree_min = float(cfg_par['rfi']['noise_measure_edges'][0])

        if (self.aperfi_rfifree_min < float(cfg_par['rfi']['lowfreq']) or self.aperfi_rfifree_min > float(cfg_par['rfi']['highfreq'])):
            self.aperfi_rfifree_min = float(cfg_par['rfi']['lowfreq'])
            self.logger.warning('\t ### Minimum frequency for RFI free range outside of observed interval ')
            self.logger.warning('\t     Correct noise_measure_edges in rfi of parameter file ###')

        self.aperfi_rfifree_max = float(cfg_par['rfi']['noise_measure_edges'][1])
        if (self.aperfi_rfifree_max < float(cfg_par['rfi']['lowfreq']) or self.aperfi_rfifree_max > float(cfg_par['rfi']['highfreq'])):
            self.aperfi_rfifree_max = float(cfg_par['rfi']['highfreq'])
            self.logger.warning('\t ### Maximum frequency for RFI free range outside of observed interval ')
            self.logger.warning('\t     Correct noise_measure_edges in rfi of parameter file ###')


        self.aperfi_rmsclip = float(cfg_par['rfi']['rms_clip'])


        chan_min = np.argmin(np.abs(self.channelFreqs[0] - self.aperfi_rfifree_min))
        chan_max = np.argmin(np.abs(self.channelFreqs[0] - self.aperfi_rfifree_max))

        time_ax_len = int(self.datacube.shape[2])

        for i in xrange(0,self.datacube.shape[0]):
            tmp_rms = np.nanmedian(self.datacube[i, chan_min:chan_max, 0])
            med2 = abs(self.datacube[i, chan_min:chan_max, 0] - tmp_rms)
            madfm = np.ma.median(med2) / 0.6744888
            flag_lim = self.aperfi_rmsclip*madfm  
            self.flag_lim_array[i] = flag_lim    

            for j in xrange(0,self.datacube.shape[1]):
                tmpar = self.datacube[i,j,:]
                mean  = np.nanmean(tmpar)
                tmpar = tmpar-mean
                tmpar = abs(tmpar)
                self.mean_array[i,j] = mean
                #change masked values to very high number
                #inds = np.where(np.isnan(tmpar))
                tmpar[np.isnan(tmpar)]=np.inf
                tmpar.sort()
                index_rms = np.argmin(np.abs(tmpar - flag_lim))
                tmp_over = len(tmpar[index_rms:-1])+1
                if tmp_over == 1. :
                    tmp_over = 0.
                rms[i,j] = 100.*tmp_over/time_ax_len

        self.write_freq_base(cfg_par,rms,time_step)

        self.logger.info('\t ... RFI found ... \n\n')


        return 0

    def write_freq_base(self,cfg_par,rms,time_step=-1) :
        '''
        
        Writes an image.fits of visibilities ordered by frequency, baseline.
        Baselines are ordered by their length.
        
        '''
        self.logger.info('\t ... Wrtiting RFI by frequency and baseline ...')


        #reverse frequencies if going from high-to-low         
        if time_step != -1:
            time_tmp = int(float(cfg_par['time_chunks']['time_step'])*time_step)
            if time_tmp == 0:
                time_name = '00'+str(time_tmp)+'m'
            elif time_tmp <100:
                time_name = '0'+str(time_tmp)+'m'
            else:
                time_name= str(time_tmp)+'m'    
        else:
            time_name = 'full' 
        
        self.rfi_freq_base = cfg_par['general']['rfidir']+'freq_base_'+time_name+'_im.fits'
        #set fits file
        hdu = fits.PrimaryHDU(rms)
        hdulist = fits.HDUList([hdu])
        header = hdulist[0].header
        #write header keywords               
        header['CRPIX1'] = 1
        header['CDELT1'] = np.abs(self.channelWidths[0,0]/1e6)
        header['CRVAL1'] = self.channelFreqs[0,0]/1e6
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
        fits.writeto(self.rfi_freq_base,rms,header,overwrite=True)
        hdulist.close()

        self.logger.info('\t ... RFI on fits file ...\n')


        return 0

    def rfi_flag(self,cfg_par):
        '''
        Creates a new MS file where RFI has been flagged in the FLAG column
        '''

        #self.flag_lim_array = np.zeros([self.vis.shape[0]])

        self.logger.info('\t ... Flagging RFI ... \n')


        for i in xrange(0,self.vis.shape[0]):
            
            if any(x == self.ant1[i] for x in self.aperfi_badant) or any(x == self.ant2[i] for x in self.aperfi_badant):
                continue
            else:
                indice=self.blMatrix[self.ant1[i],self.ant2[i]]
                flag_clip = self.flag_lim_array[indice]
            
            for j in xrange(0,self.vis.shape[1]):
                mean = self.mean_array[indice,i]
                if (np.abs(self.vis[i,j,0]) - mean ) > flag_clip:
                
                    self.flag[i,j,0] = True
        

        self.rfifile = cfg_par['general']['rfidir']+'rfi_flagged_vis.MS'

        t=tables.table(self.msfile)
        tout = t.copy(self.rfifile)       # make deep copy
        tout.close()
        t.close()
        tout = tables.table(self.rfifile, readonly=False)
        #coldmi = tout.getdminfo('DATA')     # get dminfo of existing column
        #coldmi["NAME"] = 'CORRECTED_DATA'               # give it a unique name
        #tout.addcols(tables.maketabdesc(tables.makearrcoldesc("CORRECTED_DATA",0., ndim=2)),
        #           coldmi)
        #tout.putcol('CORRECTED_DATA',vis)
        tout.putcol('FLAG',self.flag)
        tout.close()


        self.logger.info('\t ... RFI flagged ... \n\n')


        return 0

