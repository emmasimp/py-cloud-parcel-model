#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:30:36 2017

@author: mbexkes3
"""

# =============================================================================
# This file contains values of variables that will be input from a NAMELIST 
# =============================================================================

ncomps=1
nbins=70
nmodes=3
nmoms=2 # number of moments, double moment mass and number
NAER = [2400e6, 7e6, 9200e6]
D_AER = [350e-9, 1500e-9, 70e-9]
molw_aer = [556.57e-3, 556.57e-3, 132.14e-3]
sig=[0.4, 0.4, 0.43]
k=[0.0061, 0.0061, 0.61]
rhoa = [2560, 2560, 1770]
RH = 0.9077
RHint = RH
T=273-25
Tint = T
P=99500
Pint = P
Dlow=10e-9
rkm=2
dt = 1
w = 5

#constants for chamber T and P fits
PRESS1 = 6.2812e2
PRESS2 = 2.6e-3
Twall = T
Therm_chamber = 0.0015
Temp1 = 7.5
Temp2 = 0.012