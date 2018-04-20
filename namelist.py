#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:30:36 2017

@author: mbexkes3
"""

# =============================================================================
# This file contains values of variables that will be input from a NAMELIST 
# =============================================================================

nbins	=70
nmodes	=2
nmoms	=2 # number of moments, double moment mass and number, move this to variables.py
NAER 	=  [1000e6, 	1e6, 		9200e6, 	0, 0]
D_AER 	=  [100e-9, 1000e-9,	70e-9, 	0, 0]
sig	    =  [0.74,	0.65,		0.43, 		0, 0]

Mass_frac = {'ammonium sulphate': [1.0, 0.01],
            'sea salt':           [0.0,   0.0],
            'sulphuric acid':     [0,   0.0],
            'fulvic acid':        [0,   0],
            'Kaolinite':          [0,   0.0],
            'Montmorinillite':    [0.0,   0.0],
            'Feldspar':           [0.0,   0.99],
            'Illite' :            [0,   0.0],
            'Bio' :               [0,   0.0],
            'test' :              [0, 0]} # to add another aerosol type change 'test' 
                                          # to name of new aerosol type

RH = 0.95
T = 250
P = 99500
w = 10
runtime = 25

# choose type of simulation
Simulation_type = 'parcel' # parcel or chamber

# choose criteria for heteogeneous freezing
Heterogeneous_freezing_criteria = 'RH' # activation or RH

Dlow=10e-9
rkm=1
dt = 1

#constants for chamber T and P fits
PRESS1 = 6.2812e2
PRESS2 = 2.6e-3
Twall = T
Therm_chamber = 0.0015
Temp1 = 6.9e0
Temp2 = 0.03
