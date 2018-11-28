#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:56:34 2018

@author: mbexkes3
"""
from numpy import zeros

import namelist 
import constants

k = zeros(nmodes)
rhoa = zeros(nmodes)
molw_aer = zeros(nmodes)

for mode in range(nmodes):
    k[mode] = sum([Mass_frac[key][mode]*aerosol_dict[key][3] for key in aerosol_dict.keys()])
    rhoa[mode] = sum([Mass_frac[key][mode]*aerosol_dict[key][1] for key in aerosol_dict.keys()])
    molw_aer[mode] = sum([Mass_frac[key][mode]*aerosol_dict[key][0] for key in aerosol_dict.keys()])
