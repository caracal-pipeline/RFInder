#!/usr/bin/env python

import os
import logging

    def make_psf(cfg_par) :
        '''
        use wsclean to predict the psf of the observation
        '''

        command1 = '''wsclean -name psfonly -mem 100 -no-dirty -weight natural -super-weight 1.0'''
        command2 = '''-weighting-rank-filter-size 16 -size 512 512 -scale 3.0asec -channels-out 1'''
        command3 = ''' -grid-mode kb -kernel-size 7 -oversampling 63 -make-psf-only -pol xx -intervals-out 1'''
        command4 = ''' -data-column DATA -gain 0.1 -mgain 1.0 -multiscale-scale-bias 0.6 -fit-beam -no-fit-beam '''+cfg_par.rfifile

        command = command1+command2+command3+command4

        os.system(command)

        logging.info("\tPSF found\t")

        return 0