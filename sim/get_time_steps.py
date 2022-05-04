#!/usr/bin/python
#get_time_steps.py

# flake8: noqa

import pandas as pd
import os

def get_time_steps():
    # get timesteps
    if os.getcwd().split('/')[1] == 'home':
        steps = pd.read_csv(('/home/deth/Desktop/CAL/research/projects/sim/'
                             'ch2/climate_change_adaptation_and_genomic_arch/sim/'
                             'time_steps.csv'))
    # or else get filepaths on Savio
    else:
        steps = pd.read_csv(('/global/scratch/users/drewhart/'
                             'ch2/climate_change_adaptation_and_genomic_arch/sim/'
                             'time_steps.csv'))
    # set time when environmental change begins
    change_T = int(steps[steps['name']=='start']['num'].values[0])
    # set time when environmental change ends
    T = int(steps[steps['name']=='end']['num'].values[0])

    return change_T, T
