#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:30:36 2017

@author: mbexkes3
"""

# =============================================================================
# This file contains values of variables that will be input from a NAMELIST 
# =============================================================================

nbins	= 70
nmodes	= 2
ncomps = 9
n_sv = 1
#nmoms	= 2 # number of moments, double moment mass and number, move this to variables.py, not sure this variable does anything any more ...
NAER 	=  [900e6, 	55e6, 		9200e6, 	0, 0]
D_AER 	=  [350e-9, 1500e-9,	222e-9, 	0, 0]
sig	    =  [0.45,	0.4,		0.43, 		0, 0]

Mass_frac = {'ammonium sulphate': [0.5, 0.1],
            'sea salt':           [0.0,   0.0],
            'sulphuric acid':     [0,   0.0],
            'fulvic acid':        [0.0,   0.0],
            'Kaolinite':          [0,   0.0],
            'Montmorinillite':    [0.0,   0.0],
            'Feldspar':           [0.0,   0.8],
            'Illite' :            [0,   0.0],
            'Bio' :               [0,   0.0],
            'test' :              [0, 0],
            'SV01':               [0.5,0.05],
            'SV02':               [0.0,0.05],
            'SV03':               [0.0,0.0]} # to add another aerosol type change 'test' 
                                          # to name of new aerosol type

sv_mass_frac = {'SV01': [0.5,0.05],#this important! used in recalc_kappa
                'SV02': [0.0,0.05],
                'SV03': [0.0,0.0]}

RH = 0.95
T = 280
P = 100700
w = 5
runtime = 4

# choose type of simulation
simulation_type = 'chamber' # parcel or chamber

# semi-volatile flag
SV_flag = False

# choose criteria for heteogeneous freezing
Heterogeneous_freezing_criteria = 'RH>1' # activation, RH or critical mass of water
alpha_crit = 70 # for critical mass of water for freezing criteria

Dlow=10e-9
rkm=1
dt = 1

#constants for chamber T and P fits
PRESS1 = 628.12
PRESS2 = 2.6e-3
Twall = T
Therm_chamber = 0.0015
Temp1 = 7.5
Temp2 = 0.01

output_filename = 'output_file1'
