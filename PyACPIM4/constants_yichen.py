#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  8 16:25:58 2017

@author: mbexkes3
"""

# =============================================================================
# This file contains all the constants required for pyACPIM calculations 
# =============================================================================

rhow=1000
rhoi = 910
mw = 18e-3
R = 8.3144521 #joules mol-1 K-1 the universal gas constant
surface_ten=0.72e-3
JOULES_IN_A_CAL = 4.187
JOULES_IN_AN_ERG = 1e-7
ALPHA_COND = 1e0
ALPHA_THERM = 1e0
ALPHA_DEP = 1e0
ALPHA_THERM_ICE = 1e0
LV = 2.5e6
LS = 2.837e6
g = 9.8
w=5
RA = 287.05 # joules kg-1 K-1 air constant for dry air   P = rhoRT where R is the individual gas constant = Ru/Mgas
RV = 461.51 # gas constant water vapour joules kg-1 K-1
CP = 1005
CPV = 1870
CPW = 4.27e3
CPI = 2104.6e0
eps = RA/RV

# heteorgeneous ice nucleation values
ns_dict = {'Kaolinite':[0, -0.8881, 12.9445, -0.8802, 231.383],
           'Montmorinillite':[0, -0.517, 8.934, 0, -1],
           'Feldspar':[2, -0.1963, 60.2118],
           'Illite':[3, 6.53043e4, -8.2153088e2, 3.446885376, -4.82268],
           'Bio':[0, -0.89441, 19.0736],
           'None':0,
           'Contents':['equation type','ns constant 1','ns constant 2',
                       'ns constant 3','ns constant 4']}

# aerosol types
aerosol_dict = {'ammonium sulphate': [132.14e-3, 1770, 3, 0.61],
                'ammonium nitrate': [80.044e-3, 1725, 2, 0.67],
				'test': [1770,2560,2,0.6],
				'Black carbon': [200e-3, 1200, 1, 1e-5],             
                'sea salt': [58.44e-3, 2160, 2, 1.28],
                'sulphuric acid': [98.08e-3, 1840, 2.5, 0.9],
                'fulvic acid': [308.244e-3, 1500, 1, 0.067],
                'Kaolinite': [258.16e-3, 2650, 0, 0],
                'Montmorinillite': [549.07e-3, 2080, 0, 0],
                'Feldspar': [556.57e-3, 2560, 0, 0.0],
                'Illite' : [389.34e-3, 2770, 0, 0],
                'Bio' : [4700e-3, 1200, 0, 0],
                'SVC20': [200e-3, 1500, 1, 0.2, -6, 150e3],
                'SVC1': [200e-3, 1500, 1, 0.18, -5, 150e3],
                'SVC100': [200e-3, 1500, 1, 0.16, -4, 150e3],
                'SVC400' : [200e-3, 1500, 1, 0.12, -3, 150e3],
                'SV02': [190e-3, 1500, 1, 0.12, -2, 150e3],
                'SV03': [180e-3, 1500, 1, 0.1, -1, 150e3],
                'SV04': [170e-3, 1500, 1, 0.1, 0, 150e3],
                'SV05': [160e-3, 1500, 1, 0.08, 1, 150e3],
                'SV06': [150e-3, 1500, 1, 0.08, 2, 150e3], #[0.2, 0.18, 0.16, 0.14, 0.12, 0.10, 0.08, 0.06, 0.04, 0.02]
                'SV07': [140e-3, 1500, 1, 0.06, 3, 150e3]} #[1280, 1360, 1300, 1534, 1500, 1500, 1100, 1150, 1231, 1100]
                #'Contents' :['molecular weight (kg/mol)','density (kg/m3)'
                        #   nu  ,kappa, 'CSTAR',delta h vap]}
                
                #'SVC20':   [160.17e-3,   1280, 1, 0, 1.2304, 150e3],
                #'SVC1':    [146e-3,      1360, 1, 0, -0.0716, 150e3],
                #'SVC100':  [132.1146e-3, 1300, 1, 0, 2.0253,  150e3],
                #'SVC400' : [144.13e-3,   1534, 1, 0, 2.5821,  150e3]
