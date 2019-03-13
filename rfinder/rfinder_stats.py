#!/usr/bin/env python

# Import modules
import sys, string, os
import numpy as np
import yaml
import json
import glob
import casacore.tables as tables
import logging


from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn


import rfinder


class rfi_stats:



    def __init__(self):
        
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


    def make_psf(self,cfg_par) :
        '''
        use wsclean to predict the psf of the observation
        '''
        command1 = '''wsclean -name psfonly -mem 100 -no-dirty -weight natural -super-weight 1.0'''
        command2 = '''-weighting-rank-filter-size 16 -size 512 512 -scale 3.0asec -channels-out 1'''
        command3 = ''' -grid-mode kb -kernel-size 7 -oversampling 63 -make-psf-only -pol xx -intervals-out 1'''
        command4 = ''' -data-column DATA -gain 0.1 -mgain 1.0 -multiscale-scale-bias 0.6 -fit-beam -no-fit-beam '''+cfg_par['general']['rfidir']

        command = command1+command2+command3+command4

        os.system(command)

        logging.info("\tPSF found\t")

        return 0


    def time_chunk(self,cfg_par):

        self.logger.info("\t ...  Observing time Info ... \n")
        self.msfile = cfg_par['general']['msfullpath']
        
        t=tables.table(self.msfile)
        self.time = t.getcol('TIME')
        t.close()

        starttime= self.time[0]
        
        endtime=self.time[-1]
        time_chunk = float(cfg_par['rfi']['chunks']['time_step'])*60.
        
        times=np.arange(starttime,endtime+time_chunk*1.,time_chunk)
        
        cfg_par['rfi']['times'] = times
        startdate=Time(starttime/3600./24.,format='mjd',scale='utc')
        cfg_par['rfi']['startdate'] = startdate
        startdate.format='iso' 
        startdate.subformat='date_hm'       

        enddate=Time(endtime/3600./24.,format='mjd',scale='utc')
        cfg_par['rfi']['enddate'] = enddate

        enddate.format='iso'        
        enddate.subformat='date_hm'       


        self.logger.info('\t Start date: {0:%y}{0:%b}{0:%d}:{0:%X}'.format(startdate.datetime))
        self.logger.info('\t End date  : {0:%y}{0:%b}{0:%d}:{0:%X} \n\n'.format(enddate.datetime))
        
        return times,startdate,enddate

    def baseline_stats(self,cfg_par):

        self.logger.info("\t ... Baseline info ... \n")


        baseline_cutoff = float(cfg_par['rfi']['baseline_cut'])
        lenghts = np.array([cfg_par['rfi']['baseline_lenghts']])+0.
        index_baselines = (np.abs(lenghts - baseline_cutoff)).argmin()

        nrBaseline = cfg_par['rfi']['number_baseline']
        self.logger.info('\t\t Maximum baseline length     m = '+ str(np.round(lenghts[0][-1],0)))
        self.logger.info('\t\t Minimum baseline length     m = '+ str(np.round(lenghts[0][0],0)))
        self.logger.info('\t\t Total number of baselines     = '+ str(nrBaseline))


        num_short = index_baselines+1
        num_long = nrBaseline-num_short
        cfg_par['rfi']['num_short'] = num_short 
        cfg_par['rfi']['num_long'] = num_long

        self.logger.info('\t\t Number of baselines < {0:d} m = {1:d}'.format(cfg_par['rfi']['baseline_cut'],int(num_short)))
        self.logger.info('\t\t Number of baselines > {0:d} m = {1:d}\n'.format(cfg_par['rfi']['baseline_cut'],int(num_long)))


        #outdir = cfg_par['general']['rfidir']
        #basefile = outdir+str(cfg_par['rfi']['telescope'])+'_baselines.txt'

        #if os.path.exists(basefile) == False:
        #    ascii.write(baseline_array, overwrite=True)

        #self.logger.info("\t ... Baselines saved to file ... \n")

        return 0

    def predict_noise(self,cfg_par,channelWidths,interval,flag):
        '''
        
        Predict expected natural rms adding up all unflaged integrations and polarisations

        this module is takend from noisy.py [https://github.com/paoloserra/noisy]

        '''
        self.logger.info("\t ...  Predicting natural r.m.s. ... \n")
        self.msfile = cfg_par['general']['msfullpath']
        
        # Derive quantities
        kB=1380.6  

        STOKES = ['xx','yy','XX','YY','xy','yx','XY','YX']
        if cfg_par['rfi']['polarization'] in STOKES:                                # Boltzmann constant (Jy m^2 / K)
            polnum = 1
        else:
            polnum = 2
        cfg_par['rfi']['polnum'] = polnum

        tele= cfg_par['general']['telescope']

        if tele == 'meerkat' or tele == 'MeerKAT' or tele == 'meerKAT' or tele == 'meer':
            tsyseff = 30.             # Tsys/eff(K)
            diam = 13.5               # m
        elif tele == 'apertif' or tele == 'Apertif' or tele == 'APERTIF' or tele == 'wsrt':
            tsyseff = 93.
            diam = 25.
        else:
            self.logger.error('\t Telescope not known or not specified, will use WSRT instead')
            tsyseff = 93.
            diam = 25.

        Aant=np.pi*(diam/2)**2                      # collecting area of 1 antenna (m^2)
        SEFD=2*kB*tsyseff/Aant                  # frequency independent system equivalent flux density (Jy)
        
        # Print assumptions
        self.logger.info('\t Assumptions on {0:s} telescope'.format(tele))
        self.logger.info('\t\tDish diameter = {0:.1f} m'.format(diam))
        self.logger.info('\t\t ... and SEFD = {0:0f} Jy'.format(SEFD))
        self.logger.info('\t\t ... and Tsys = {0:.1f} K'.format(tsyseff))

        self.logger.info('\t Properties of observation'.format(tele))
        nrBaseline = cfg_par['rfi']['number_baseline']
        self.logger.info('\t\t Total number of baselines = '+ str(nrBaseline))
        self.logger.info('\t\t Total number of channels = '+ str(channelWidths.shape[1]))
        cfg_par['rfi']['total_channels'] = channelWidths.shape[1]
        cfg_par['rfi']['exptime'] = interval.sum()/nrBaseline/3600.
        self.logger.info('\t\t Observing time on source = {0:.5f} h ({1:d} polarisations)\n'.format(cfg_par['rfi']['exptime'],polnum))

        rms=np.sqrt(2)*kB*tsyseff/Aant/np.sqrt(channelWidths*interval.sum()*polnum)

        if len(rms.shape)==2 and rms.shape[0]==1: rms=rms[0]

        self.logger.info('\t Stokes I natural r.m.s.       = {0:.3e} mJy/b '.format(np.nanmedian(rms*1e3)))

        cfg_par['rfi']['theo_rms'] = rms
        
        self.logger.info("\t ... Natural r.m.s. predicted ... \n")


        return 0


    def alt_az(self,cfg_par,time):

        self.logger.info("\t ... Altitude/Azimuth info ... \n")

        #open miriad observation
        tele= cfg_par['general']['telescope']
        coord = cfg_par['rfi']['coords']

        if tele == 'meerkat' or tele == 'MeerKAT' or tele == 'meerKAT' or tele == 'meer':
            telescope = EarthLocation(lat= -30.712856*u.deg, lon=21.44377374*u.deg, height=1051.0*u.m)        
        elif tele == 'apertif' or tele == 'Apertif' or tele == 'APERTIF' or tele == 'wsrt':
            # Set position of Westerbork (source: http://www.sonel.org/spip.php?page=gps&idStation=884)
            telescope = EarthLocation(lat= 52.91460037*u.deg, lon=6.60449982*u.deg, height=82.2786*u.m)        
        else:
            self.logger.error('\t Observing location unknown, impossible to estimate Alt/Az')
        
        time = Time(time/3600./24.,format='mjd',scale='utc')
        frame = AltAz(obstime=time, location=telescope)
        obs_altaz = coord.transform_to(frame)
        cfg_par['rfi']['altaz'] = obs_altaz
       
        self.logger.info(('\t\t ... Altitude = {0:s}').format(cfg_par['rfi']['altaz'].alt))
        self.logger.info(('\t\t ... Azimuth = {0:s}\n').format(cfg_par['rfi']['altaz'].az))

        self.logger.info("\t ... Alt/Az done ... \n")

        return obs_altaz