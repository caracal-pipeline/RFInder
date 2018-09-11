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
import numpy as np
from astropy import units as u

from astropy.io import fits, ascii
from astropy import units as u
from astropy.table import Table, Column, MaskedColumn
import rfinder


class rfi_stats:



    def __init__(self):

        #set logger
        logging.basicConfig(level=logging.INFO)
        self.logger = logging.getLogger(__name__)

    def make_psf(cfg_par) :
        '''
        use wsclean to predict the psf of the observation
        '''
        print cfg_par['general']['rfidir']
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


    def predict_noise(self,cfg_par,channelWidths,interval,flag):
        '''
        
        Predict expected natural rms adding up all unflaged integrations and polarisations

        this module is takend from noisy.py [https://github.com/paoloserra/noisy]

        '''
        self.logger.info("\t ...  Predicting natural r.m.s. ... \n")
        self.msfile = cfg_par['general']['msfullpath']
        
        # Derive quantities
        kB=1380.6                                   # Boltzmann constant (Jy m^2 / K)
        tele= cfg_par['rfi']['telescope']

        if tele == 'meerkat' or tele == 'MeerKAT' or tele == 'meerKAT' or tele == 'meer':
            tsyseff = 30.             # Tsys/eff(K)
            diam = 13.5               # m
            polnum = 4
        elif tele == 'apertif' or tele == 'Apertif' or tele == 'APERTIF' or tele == 'wsrt':
            tsyseff = 93.
            diam = 25.
            polnum = 2
        else:
            self.logger.error('\t Telescope not known or not specified, will use WSRT instead')
            tsyseff = 93.
            diam = 25.
            polnum = 2

        Aant=np.pi*(diam/2)**2                      # collecting area of 1 antenna (m^2)
        SEFD=2*kB*tsyseff/Aant                  # frequency independent system equivalent flux density (Jy)
        
        # Print assumptions
        self.logger.info('\t Assumptions on {0:s} telescope'.format(tele))
        self.logger.info('\t\tDish diameter = {0:.1f} m'.format(diam))
        self.logger.info('\t\t ... and SEFD = {0:0f} Jy'.format(SEFD))
        self.logger.info('\t\t ... and Tsys = {0:.1f} K'.format(tsyseff))
        self.logger.info('\t\t antenna diam = {0:.1f} m\n'.format(diam)) 

        self.logger.info('\t Assumptions on {0:s} telescope'.format(tele))
        nrBaseline = cfg_par['rfi']['number_baseline']
        self.logger.info('\t\t Total number of channels = '+ str(channelWidths.shape[1]))
        self.logger.info('\t\t Observing time on source = {0:.5f} h ({1:d} polarisations)\n'.format(interval.sum()/nrBaseline/3600.,polnum))
        
        rms=np.sqrt(2)*kB*tsyseff/Aant/np.sqrt(channelWidths*interval.sum()*4)

        if len(rms.shape)==2 and rms.shape[0]==1: rms=rms[0]

        self.logger.info('\t Stokes I natural rms       = {0:.3e} mJy/b '.format(np.nanmedian(rms*1e3)))

        cfg_par['rfi']['theo_rms'] = rms

        self.logger.info("\t ... Natural r.m.s. predicted ... \n")


        return 0