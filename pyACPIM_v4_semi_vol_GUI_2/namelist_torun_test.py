#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:30:36 2017

@author: mbexkes3
"""

# =============================================================================
# This file contains values of variables that will be input from a NAMELIST 
# =============================================================================

nbins	= 60
nmodes = 1 #	= 2
ncomps = 9
n_sv = 1

NAER = [300e6,600e6,900e6]# 	=  [200e6, 	2.19e8, 		9200e6, 	0, 0]
D_AER = [600e-9,600e-9,200e-9]# 	=  [17e-9, 95e-9,	222e-9, 	0, 0]
sig = [0.5,0.53,0.78]#	    =  [0.3,	0.3,		0.43, 		0, 0]

Mass_frac = {'ammonium sulphate': [0.5,   0.000, 0.000],
            'sea salt':           [0.0,   0.0,0.0],
            'sulphuric acid':     [0.0,   0.0,0.0],
            'fulvic acid':        [0.0,   0.0,0.0],
            'Kaolinite':          [0.0,   0.0,0.0],
            'Montmorinillite':    [0.0,   0.0,0.0],
            'Feldspar':           [0.0,   0.0,0.0],
            'Illite' :            [0.0,   0.0,0.0],
            'Bio' :               [0.0,   0.0,0.0]}
            

semi_vols = ['SVC20','SVC100']                                      
SV_MF = [0.5,0.0,0.0]#[0.25,0.25] # this needs to be [:n_sv]
SV_MR = [1.9e-10,1.9e-10,0.053e-9]# = 1.39e-12 # this needs to be [:n_sv]
SV_flag = True
kappa_flag = True


RH = 0.99# 0.8383
T = 280.0#289.5#289.5#287.9217
P = 95000# 97693
w = 1
runtime = 2

# choose type of simulation
simulation_type = 'parcel' # parcel or chamber

# choose criteria for heteogeneous freezing
Heterogeneous_freezing_criteria = 'RH>1' # activation, RH or critical mass of water
alpha_crit = 70 # for critical mass of water for freezing criteria

Dlow=10e-9
rkm=1
dt = 1

#constants for chamber T and P fits
PRESS1 = 512.3# = 512.3
PRESS2 = 0.0044068# = 0.0044068
Twall = T
Therm_chamber = 0.0015
Temp1 = 10.484769# = 6.36909229# = 8.4179284# = 8.3
Temp2 = 0.0357931506# = 0.0401# = 0.037410977# = 0.025

output_filename = 'test.nc'# = 'output_run4_C400_T_U.nc'# = 'output_run4_C400_T_mean.nc'# = 'output_sv1_run3_2403.nc'
