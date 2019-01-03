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
nmodes	= 1
#nmoms	= 2 # number of moments, double moment mass and number, move this to variables.py, not sure this variable does anything any more ...
NAER 	=  [793688447.333, 	1100881232.88, 		9200e6, 	0, 0]
D_AER 	=  [2.90764523934e-07, 4.29312369463e-08,	222e-9, 	0, 0]
sig	    =  [0.408076557823,	0.554713826442,		0.43, 		0, 0]

Mass_frac = {'ammonium sulphate': [0.1, 0.01],
            'sea salt':           [0.0,   0.0],
            'sulphuric acid':     [0,   0.0],
            'fulvic acid':        [0.0,   0.0],
            'Kaolinite':          [0,   0.0],
            'Montmorinillite':    [0.0,   0.0],
            'Feldspar':           [0.0,   0.99],
            'Illite' :            [0,   0.0],
            'Bio' :               [0,   0.0],
            'test' :              [0, 0]} # to add another aerosol type change 'test' 
                                          # to name of new aerosol type

RH = 0.759975932311
T = 223.677609662
P = 96000
w = 1.5
runtime = 15

# choose type of simulation
Simulation_type = 'parcel' # parcel or chamber

# choose criteria for heteogeneous freezing
Heterogeneous_freezing_criteria = 'RH' # activation, RH or critical mass of water
alpha_crit = 70 # for critical mass of water for freezing criteria

Dlow=10e-9
rkm=1
dt = 1

#constants for chamber T and P fits
PRESS1 = 628.12
PRESS2 = 2.6e-3
Twall = T
Therm_chamber = 0.0015
Temp1 = 6.9
Temp2 = 0.03

output_filename = 'output_file'
