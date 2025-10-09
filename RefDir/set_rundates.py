#!/usr/bin/env python3

import os
from datetime import datetime, timedelta
import shutil
import fileinput
import sys

targ_yr = int(sys.argv[1])
start_month = int(sys.argv[2])
dir_str=os.path.basename(os.getcwd())
sim_type = dir_str.split('_')[-2]
if sim_type == 'RF':
    stop_yr = targ_yr
    stop_month = start_month
    stop_day = 2
elif sim_type == 'Std':
    stop_yr = targ_yr + 1
    stop_month = 1
    stop_day = 1
else:
    raise ValueError('Directory name not correctly configured')

with fileinput.FileInput('input.geos', inplace=True, backup='.bak') as file:
    for line in file:
        if line.startswith('Start YYYYMMDD'):
            subline = line.split(':')
            line = subline[0] + ': {:04d}{:02d}01 000000\n'.format(targ_yr,start_month)
        elif line.startswith('End   YYYYMMDD'):
            subline = line.split(':')
            line = subline[0] + ': {:04d}{:02d}{:02d} 000000\n'.format(stop_yr,stop_month,stop_day)
        print(line,end='')
