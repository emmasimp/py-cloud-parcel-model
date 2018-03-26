#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:30:36 2017

@author: mbexkes3
"""

# =============================================================================
# This file contains values of variables that will be input from a NAMELIST 
# =============================================================================

ncomps	=9
nbins	=70
nmodes	=2
nmoms	=2 # number of moments, double moment mass and number
NAER 	=  [240e6, 	7e6, 		9200e6, 	0, 0]
D_AER 	=  [350e-9, 	1500e-9,	70e-9, 		0, 0]
molw_aer = [132.14e-3, 	556.57e-3,	556.57e-3,	0, 0]
sig	=  [0.4,	0.4,		0.43, 		0, 0]
k	=  [0.61, 	0.0061, 	0.0061, 	0, 0]
rhoa = 	   [1770,	2560,		2560,     	0, 0]
ns_1 =     [0,		-0.1963,	-0.1963,	0, 0]
ns_2 =     [0,		60.2118,	60.2118,	0, 0]
RH = 0.8
RHint = RH
T=260
Tint = T
P=99500
Pint = P
Dlow=10e-9
rkm=1
dt = 1
w = 5
runtime = 200

#constants for chamber T and P fits
PRESS1 = 6.2812e2
PRESS2 = 2.6e-3
Twall = T
Therm_chamber = 0.0015
Temp1 = 6.9e0
Temp2 = 0.03
