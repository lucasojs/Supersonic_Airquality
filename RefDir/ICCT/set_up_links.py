#!/usr/bin/env python3

import os

base_dir_in = '/net/d13/data/inessanz/SST_AEIC/GEOS-Chem_Emission/'
base_dir_out = '/net/d05/data/seastham/GEOS-Chem/Supersonic_2021/ICCT_Emissions'

restriction_case = 'Unrestricted'
fuel_case = 'JetA1'
for icct_case in [1,2,3,4]:
    dir_in = os.path.join(base_dir_in,'Supersonic_icct{:d}'.format(icct_case))
    if icct_case == 1:
        aircraft_case = 'BizJet'
        range_case = 4500
    elif icct_case == 2:
        aircraft_case = 'BoomCA'
        range_case = 4500
    elif icct_case == 3:
        aircraft_case = 'BizJet'
        range_case = 4000
    elif icct_case == 4:
        aircraft_case = 'BoomCA'
        range_case = 4000
    dir_out = os.path.join(base_dir_out,aircraft_case,fuel_case,'R{:d}'.format(range_case),restriction_case)
    if os.path.isdir(dir_out):
        os.rmdir(dir_out)
    os.symlink(dir_in,dir_out,target_is_directory=True)
