#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:51:12 2017

@author: mbexkes3
"""
import numpy as np
#import scipy
from scipy.misc import derivative
import matplotlib.pyplot as plt

from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

import constants as c
import namelist as n
import functions as f
import setup_grids as s

import importlib
importlib.reload(n)
importlib.reload(f)
importlib.reload(s)

IPRESS = 0 + n.nbins*n.nmodes
ITEMP = 1 + n.nbins*n.nmodes
IRH = 2 + n.nbins*n.nmodes

Y = np.zeros(n.nbins*n.nmodes+3)

Y[IPRESS] = n.Pint
Y[ITEMP] = n.Tint
Y[IRH] = n.RHint

rhobin = np.zeros(n.nmodes*n.nbins)
kappabin = np.zeros(n.nmodes*n.nbins) # just externally mixed aerosol at the moment, need to add mass fraction term for internally mixed aerosol
molwbin = np.zeros(n.nmodes*n.nbins)

for i in range(n.nmodes):
    rhobin[i*n.nbins:n.nbins+n.nbins*i] = n.rhoa[i]
    kappabin[i*n.nbins:n.nbins+n.nbins*i] = n.k[i]
    molwbin[i*n.nbins:n.nbins+n.nbins*i] = n.molw_aer[i]
    
############################### set-up grid ###################################
GRID = s.setup_grid(rhobin, kappabin, n.rkm, n.Dlow, n.nbins, n.nmodes, n.RH, n.sig, n.NAER, n.D_AER, n.T)

DIAM_EDGES = GRID[3] # does not change
MWATGRID = GRID[1] # does not change
mwat_centre =GRID[2] # this changes after rebinning
n_aer_bin = GRID[0] # the changes after rebinning
AER_MASS_CENTRE = GRID[5] # does not change as no semi-vols or coagulation
AER_MASS_EDGES = GRID[4] # does not change as no semi-vols or coagulation

Y[0:n.nbins*n.nmodes] = mwat_centre

###############################################################################

    
def dy_dt_func(t,Y):
    
    dy_dt = np.zeros(len(Y))
    #pressure, hydrostatic equation
    dy_dt[IPRESS] = -(Y[IPRESS]*c.g*n.w)/(c.RA*Y[ITEMP])
    
    #temperature change due to expansion and condensation
    #no freezing
    WV = c.eps*Y[IRH]*f.svp_liq(Y[ITEMP])/(Y[IPRESS] - f.svp_liq(Y[ITEMP]))
    WL = sum(n_aer_bin*Y[0:n.nbins*n.nmodes])
    Rm = c.RA+WV*c.RV
    Cpm = c.CP + WV*c.CPV + WL*c.CPW # just for vapour and liquid, no ice yet
    dy_dt[ITEMP] = (Rm*Y[ITEMP]*dy_dt[IPRESS])/(Cpm*Y[IPRESS])- c.LV*WV/Cpm
    
    #change in vapour content:
    # 1. equilibruim size of particles
    KK01 = f.kk01(Y[0:n.nbins*n.nmodes], Y[ITEMP], AER_MASS_CENTRE, rhobin, kappabin, molwbin)
    Dw = KK01[2] # wet diameter
    RHOAT = KK01[1] # density of particles inc water and aerosol mass
    RH_EQ = KK01[0] # equilibruim diameter
    Dd = KK01[3]
    
    # 2. growth rate of particles, Jacobson p455
    growth_rate = f.DROPGROWTHRATE(Y[ITEMP],Y[IPRESS],Y[IRH],RH_EQ,RHOAT,Dw) # rate of change of radius
    
    # 3. Mass of water condensing 
    dy_dt[0:n.nbins*n.nmodes] = np.pi*RHOAT*Dw**2 * growth_rate # change in mass of water per particle
   
    # 4. Change in vapour content
    dwv_dt = -1*sum(n_aer_bin*dy_dt[0:n.nbins*n.nmodes]) # change in water vapour mixing ratio
    
    #RH change
    dy_dt[IRH] = f.svp_liq(Y[ITEMP])*dwv_dt*(Y[IPRESS]-f.svp_liq(Y[ITEMP]))
    dy_dt[IRH] = dy_dt[IRH] + f.svp_liq(Y[ITEMP])*WV*dy_dt[IPRESS]
    dy_dt[IRH] = dy_dt[IRH] - WV*Y[IPRESS]*derivative(f.svp_liq,Y[ITEMP],dx=1.0)*dy_dt[ITEMP]
    dy_dt[IRH] = dy_dt[IRH]/(c.eps*f.svp_liq(Y[ITEMP])**2)
    
    dy_dt_func.wv = WV
    dy_dt_func.wl = WL
    dy_dt_func.rt = WV+WL
    dy_dt_func.rheq = RH_EQ
    dy_dt_func.dw = Dw
     
    return dy_dt

#--------------------- SET-UP solver ------------------------------------------

y0 = Y
t0 = 0.0

#test = dy_dt_func(1,Y)

#define assimulo problem
exp_mod = Explicit_Problem(dy_dt_func,y0,t0)

# define an explicit solver
exp_sim = CVode(exp_mod)

#set parameters
tol_list = np.zeros_like(Y)
tol_list[0:n.nbins*n.nmodes] = 1e-25
tol_list[IPRESS] = 10
tol_list[ITEMP] =1e-4
tol_list[IRH] = 1e-8 # set tolerance for each dydt function
exp_sim.atol = tol_list
exp_sim.rtol = 1.0e-6
exp_sim.inith = 1.0e-2 # initial time step-size
exp_sim.usejac = False
exp_sim.report_continuously = True
exp_sim.maxncf = 1000 # max number of convergence failures allowed by solver


#run simulation
t_final = 200
outputs = range(0, t_final, 1)
t_output,y_output = exp_sim.simulate(t_final,ncp_list=outputs)

############################# testing #########################################
plt.plot(y_output[:,IRH])
plt.show()
