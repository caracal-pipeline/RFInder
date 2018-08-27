import sys, string, os
import numpy as np

import rfinder as rfinder

file = sys.argv[1]
cfg = open(file)

rfi_par=rfinder.rfinder(file)

run = rfi_par.go(rfi_par.cfg_par)

if run == 0: 
    print '\t+------+\n\t Done \n\t+------+'