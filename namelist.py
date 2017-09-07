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
nbins=60
nmodes=2
nmoms=2 # number of moments, double moment mass and number
NAER = [2e6, 100e6]
D_AER = [200e-9, 2000e-9]
molw_aer = [132.14e-3, 132.14e-3]
sig=[0.5, 0.8]
k=[0.61, 0.61]
rhoa = [1770,1770]
RH = 0.9
RHint = RH
T=273+30
Tint = T
P=100000
Pint = P
Dlow=10e-9
rkm=1
w = 5

#constants for chamber T and P fits
PRESS1 = 6.038e2
PRESS2 = 2.947e-3
Twall = T
Therm_chamber = 0.0075
Temp1 = 6.233
Temp2 = 8.33e-3