#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  4 14:54:15 2018

this model will create, write and save output from pyACPIM
as a netcdf file

@author: emma
"""

from netCDF4 import Dataset
import namelist as n
import numpy as np

class output_nc(object):
    
    def __init__(self):
        self = Dataset(n.output_filename+'.nc','w', format='NETCDF4_CLASSIC')
        
        #create dimensions
        self.createDimension('times',n.runtime/n.dt)
        self.createDimension('bins',n.nbins)
        
        #create variables
        self.createVariable('total_ice',np.float32, ('times','bins'))
        
        
    def write_to_netcdf(self,idx, output_ice_test):
        
       self['total_ice'][idx,:] = np.sum(output_ice_test[:]) 
       
     