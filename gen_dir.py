#!/usr/bin/env python3
import os
import sys
from datetime import datetime, timedelta

if len(sys.argv) != 5:
    raise ValueError('Need exactly 4 arguments: CH4 type, scenario name, simulation round, and whether to use WACCM aerosol')

ch4_type  = sys.argv[1]
scen_name = sys.argv[2]
sim_round = int(sys.argv[3])
use_WACCM = sys.argv[4] == 'T'

if sys.argv[4] not in ['T','F']:
    raise ValueError('WACCM option must be T or F')

assert ch4_type in ['Fixed','Flux'], 'Invalid CH4 type'
is_baseline = 'Baseline' in scen_name

if use_WACCM:
    WACCM_str = ''
else:
    WACCM_str = 'noWACCM'

base_dir = os.getcwd()
for sim_type in ['Std','RF']:
    os.chdir(base_dir)
    is_RRTMG = sim_type == 'RF'
    dir_name  = 'Run_{}_{:s}_R{:d}'.format(scen_name,sim_type,sim_round)
    dir_inter = os.path.join('CH4_' + ch4_type,dir_name)
    dir_full = os.path.join(base_dir,dir_inter)
    assert not os.path.isdir(dir_full), 'Directory {} already exists'.format(dir_inter)
    rc = os.system('cp -a RefDir ' + dir_full)
    assert rc == 0, 'Copy failed'
    run_opts = '{:d} {:s} {} {} {}'.format(sim_round,ch4_type,is_baseline,use_WACCM,is_RRTMG)
    if is_RRTMG:
        run_opts = '{} {:d} {}'.format(run_opts,1,not is_baseline)
    os.chdir(dir_full)
    ems_dir = 'CH4Emissions_R{:d}'.format(sim_round)
    dst_dir = os.path.join(dir_full,ems_dir)
    if ch4_type == 'Flux':
        # Generate symlink to fixed dir
        src_dir = os.path.join(base_dir,'CH4_Fixed','Run_Baseline{:s}_Std_R{:d}'.format(WACCM_str,sim_round),ems_dir)
        os.symlink(src_dir,dst_dir,target_is_directory=True)
    else:
        os.mkdir(dst_dir)
    # Sort out restarts
    if is_RRTMG:
        os.remove(os.path.join(dir_full,'run_gcc.sh'))
        os.remove(os.path.join(dir_full,'batch_gcc.sh'))
        if is_baseline:
            os.mkdir(os.path.join(dir_full,'DynHeating'))
        else:
            src_dir = os.path.join(base_dir,'CH4_' + ch4_type,'Run_Baseline{:s}_RF_R{:d}'.format(WACCM_str,sim_round),'DynHeating')
            os.symlink(src_dir,os.path.join(dir_full,'DynHeating'),target_is_directory=True)
        # Link to non-RRTMG files
        std_dir = os.path.join(base_dir,'CH4_' + ch4_type,'Run_{}_{:s}_R{:d}'.format(scen_name,'Std',sim_round))
        os.symlink(os.path.join(std_dir,'restarts'),
                   os.path.join(dir_full,'restarts'))
        # Link HEMCO restart files
        for l_yr in range(2000,2014):
            for l_mo in range(1,13):
                f = 'HEMCO_restart.{:04d}{:02d}010000.nc'.format(l_yr,l_mo)
                f_src = os.path.join(std_dir,f)
                f_dst = os.path.join(dir_full,f)
                os.symlink(f_src,f_dst) 
    else:
        os.remove(os.path.join(dir_full,'run_rrtmg.sh'))
        os.mkdir(os.path.join(dir_full,'restarts'))
        os.mkdir(os.path.join(dir_full,'restarts','tracers'))
        os.mkdir(os.path.join(dir_full,'restarts','cspec'))
        if sim_round == 1:
            # Link to initial restart
            os.symlink(os.path.join(base_dir,'RefData','initial_trac_rst_2x25.merra_L72.T93'),
                       os.path.join(dir_full,'restarts','tracers','trac_rst.merra_L72.T93.200001010000'))
            #os.symlink(os.path.join(base_dir,'RefData','initial_HEMCO_rst.nc'),
            #           os.path.join(dir_full,'HEMCO_restart.200001010000.nc'))
        else:
            src_dir = os.path.join(base_dir,'CH4_' + ch4_type,'Run_{}_{:s}_R{:d}'.format(scen_name,sim_type,sim_round-1))
            os.symlink(os.path.join(src_dir, 'restarts','tracers','trac_rst.merra_L72.T93.201401010000'),
                       os.path.join(dir_full,'restarts','tracers','trac_rst.merra_L72.T93.200001010000'))
            os.symlink(os.path.join(src_dir, 'restarts','cspec','spec_rst.merra_4x5_UCX.2014010100'),
                       os.path.join(dir_full,'restarts','cspec','spec_rst.merra_4x5_UCX.2000010100'))
            os.symlink(os.path.join(src_dir, 'HEMCO_restart.201401010000.nc'),
                       os.path.join(dir_full,'HEMCO_restart.200001010000.nc'))
    rc = os.system('./setup_dir.py 2000 {}'.format(run_opts))
    assert rc == 0, 'Setup failed in ' + dir_inter + ' with run options {}'.format(run_opts)
    
