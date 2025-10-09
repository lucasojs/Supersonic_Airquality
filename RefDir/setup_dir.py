#!/usr/bin/env python3

import os
from datetime import datetime, timedelta
import shutil
import fileinput
import sys

run_opts_file = 'run_opts.dat'
opt_list = sys.argv
if len(sys.argv) == 2:
    assert os.path.isfile(run_opts_file), 'If only year is supplied, need run_opts.dat to exist'
    with open(run_opts_file,'r') as f:
        content = f.readlines()
    opt_list += [x.strip() for x in content]
elif len(sys.argv) < 7 or (sys.argv[6] == 'False' and len(sys.argv) != 7) or (len(sys.argv) != 9 and sys.argv[6] == 'True'):
    scr_name = sys.argv[0]
    #print('N', len(sys.argv))
    print('Need 5-7 options to be passed: year, round, CH4 type, reference run, WACCM, and RRTMG are needed at minimum')
    print('If using RRTMG, also need the month and if dynamic heating is being read in')
    print('Example: {} 2000 1 Fixed False True False     # First year, 1st round, standard, non-reference run with WACCM, without RRTMG and fixed CH4'.format(scr_name))
    print('         {} 2002 2 Flux True True True 3 True # Third year, 2nd round, reference, with WACCM, RRTMG sim for March with FDH read-in and flux CH4'.format(scr_name))
    print('Alternatively, if only a year is passed, script will attempt to read options from the file run_opts.dat')
    raise ValueError('Invalid inputs')

for i_opt in [4,5,6,8]:
    if len(opt_list) >= i_opt + 1:
        assert opt_list[i_opt] in ['True','False'], 'Invalid choice for option {:d}'.format(i_opt)

# Which year are we running?
targ_yr = int(opt_list[1])

# Which 14-year round is this?
sim_round = int(opt_list[2])

# Is this with fixed or flux CH4?
ch4_type = opt_list[3]

# Is this a reference simulation?
reference_run = opt_list[4] == 'True'

# Are we using WACCM aerosols and SNAP?
use_WACCM = opt_list[5] == 'True'

# Are we running a one-day RRTMG simulation?
run_rrtmg = opt_list[6] == 'True'

# Are we generating dynamic heating data or reading it in?
if run_rrtmg:
    start_month = int(opt_list[7])
    stop_month = start_month
    stop_day = 2
    stop_yr = targ_yr
    read_FDH = opt_list[8] == 'True'
    if not read_FDH:
        assert reference_run, 'DH data should only be written during a reference run. Turn off FDH read-in'
else:
    start_month = 1
    stop_month = 1
    stop_yr = targ_yr + 1
    stop_day = 1
    read_FDH = False

assert targ_yr >= 2000 and targ_yr <= 2013, 'Target year must be between 2000 and 2013 inclusive'
assert sim_round > 0, 'Simulation round must be at least 1'
assert ch4_type in ['Fixed','Flux'], 'Only Fixed and Flux are valid simulation types'
assert start_month > 0 and start_month < 13, 'Start month not valid'

first_yr = targ_yr == 2000 and sim_round == 1
start_date = datetime(targ_yr,start_month,1,0,0,0)

print('Setting up directory with the following options:')
print(' => Start date      : {:s}'.format(start_date.strftime('%Y-%m-%d %H:%M')))
print(' => Simulation round: {:d}'.format(sim_round))
print(' => CH4 type        : {:s}'.format(ch4_type))
print(' => Reference run   : {}  '.format(reference_run))
print(' => WACCM + SNAP    : {}  '.format(use_WACCM))
print(' => Use RRTMG       : {}  '.format(run_rrtmg))
if run_rrtmg:
    if read_FDH:
        FDH_action = 'read'
    else:
        FDH_action = 'write'
    print(' => Dyn heating data: {:s}'.format(FDH_action))

# Check to see that we have all the inputs we need
f_list = []
if read_FDH:
    # Need all the dynamical heating inputs
    #print('CHECKING FDH')
    doy = (datetime(targ_yr,start_month,1) - datetime(targ_yr,1,1)).days + 1
    for t in range(0,24,3):
        f_list.append(os.path.join('DynHeating','DynHeating.{:04d}{:03d}{:04d}.nc'.format(targ_yr,doy,130 + 100*t)))
if not first_yr:
    # Check for CSPEC file
    f_list.append(os.path.join('restarts','cspec','spec_rst.merra_4x5_UCX.{:s}'.format(start_date.strftime('%Y%m%d%H'))))
# Check that the tracer restart file exists
f_list.append(os.path.join('restarts','tracers','trac_rst.merra_L72.T93.{:s}'.format(start_date.strftime('%Y%m%d%H%M'))))

for f in f_list:
    if 'DynHeating' in f:
        if not (os.path.isfile(f) or os.path.islink(f)):
            print('WARNING: File {} not found'.format(f))
    else:
        assert os.path.isfile(f) or os.path.islink(f), 'Could not find {}'.format(f)

if run_rrtmg:
    # Make sure that the HEMCO configurations are identical
    dir_name = os.path.basename(os.getcwd())
    std_name = dir_name.replace('_RF_','_Std_')
    rc = os.system('diff HEMCO_Config.rc ../{:s}/.'.format(std_name))
    assert rc == 0, 'HEMCO_Config differs from standard simulation!'

if ch4_type == 'Flux':
    # Need defined emissions
    CH4_SBC = False
    CH4_read_str = 'READ'
    assert os.path.islink('CH4Emissions_R{:d}'.format(sim_round))
elif ch4_type == 'Fixed':
    CH4_SBC = True
    if reference_run:
        CH4_read_str = 'WRITE'
        assert os.path.isdir('CH4Emissions_R{:d}'.format(sim_round))
    else:
        CH4_read_str = 'OFF'

def bool_str(var):
    if var:
        return 'T'
    else:
        return 'F'

replace_opts = {'Turn on RRTMG'   : bool_str(run_rrtmg),
                'Calculate LW'    : bool_str(run_rrtmg),
                'Calculate SW'    : bool_str(run_rrtmg),
                'Clear-sky flux'  : bool_str(run_rrtmg),
                'All-sky flux'    : bool_str(run_rrtmg),
                'Apply strat adj' : bool_str(run_rrtmg),
                'Read FDH'        : bool_str(read_FDH),
                'Make new restart': bool_str(not run_rrtmg),
                'save CSPEC_FULL' : bool_str(not run_rrtmg),
                'CH4 read/write'  : CH4_read_str,
                'Use WACCM'       : bool_str(use_WACCM),
                'Turn on ND49 diagnostic': bool_str(not run_rrtmg),
                'CH4 IO folder'   : 'CH4Emissions_R{:d}'.format(sim_round)}

for mo in range(1,13):
    mo_start = datetime(2001,mo,1,0,0,0)
    if mo == 12:
        mo_end = datetime(2002,1,1,0,0,0)
    else:
        mo_end = datetime(2001,mo+1,1,0,0,0)
    n_days = (mo_end-mo_start).days
    if mo == 2:
        n_days = 29 # input.geos is set up to handle leap years
    if run_rrtmg:
        output_str = '3' * n_days
    else:
        output_str = '3' + ('0' * (n_days - 1))
    output_str = ''.join(output_str)
    mo_str = 'Schedule output for {:s}'.format(mo_start.strftime('%b').upper())
    replace_opts[mo_str] = output_str

if run_rrtmg:
    replace_opts['Radiative output'] = 'Radiative output  : 72   all'
else:
    replace_opts['Radiative output'] = 'Radiative output  :  0   all'

false_MRs = ['strat Bry',' Br species']
in_glob_MRs = False
with fileinput.FileInput('input.geos', inplace=True, backup='.bak') as file:
    for line in file:
        if line.startswith('Start YYYYMMDD'):
            subline = line.split(':')
            line = subline[0] + ': {:04d}{:02d}01 000000\n'.format(targ_yr,start_month)
        elif line.startswith('End   YYYYMMDD'):
            subline = line.split(':')
            line = subline[0] + ': {:04d}{:02d}{:02d} 000000\n'.format(stop_yr,stop_month,stop_day)
        if not in_glob_MRs:
            in_glob_MRs = 'Set initial glob MRs' in line
            # Can only check for CH4 here - don't want to overwrite global VMR switch
            if ' => CH4?' in line:
                subline = line.split(':')
                line = subline[0] + ': ' + bool_str(CH4_SBC) + '\n'
        else:
            in_glob_MRs = 'Use RCP or WMO' not in line
            if in_glob_MRs:
                subline = line.split(':')
                replace_MR = first_yr
                for spc in false_MRs:
                    replace_MR = replace_MR and spc not in subline[0]
                if replace_MR:
                    line = subline[0] + ': T\n'
                else:
                    line = subline[0] + ': F\n'
        # Look for RRTMG options
        for seeker, opt in replace_opts.items():
            if seeker in line:
                subline = line.split(':')
                line = subline[0] + ': {}\n'.format(opt)
        print(line,end='')

# Write options (except year) to file
with open(run_opts_file,'w') as f:
    f.writelines([x + '\n' for x in opt_list[2:]])
