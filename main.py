#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:51:12 2017

@author: mbexkes3
"""
import numpy as np
from scipy.misc import derivative
from scipy.optimize import minimize, minimize_scalar
from scipy.optimize import brent

import matplotlib.pyplot as plt

from assimulo.solvers import CVode
from assimulo.problem import Explicit_Problem

import constants as c
import namelist as n
import functions as f
import setup_grids_cumulative_lognormal as s

import importlib # NOTE: need to be in the correct directory for this to work
importlib.reload(n)
importlib.reload(f)
importlib.reload(s)

############# SET-UP INDEXS ###################################################
IND1 = n.nbins*n.nmodes # index for mass
IND2 = IND1*2 # index for number
IPRESS = -3
ITEMP = -2
IRH = -1

IPRESS_ICE = -3
ITEMP_ICE = -2
IRH_ICE = -1
###############################################################################
t_final = n.runtime
Y = np.zeros(n.nbins*n.nmodes*2+3)
Y_AER = np.zeros([int(t_final/n.dt),IND1]) # mass of aerosol bin
YICE_AER = np.zeros_like(Y_AER) # mass of aerosol in ice bins

Y[IPRESS] = n.Pint
Y[ITEMP] = n.Tint
Y[IRH] = n.RHint

rhobin = np.zeros(n.nmodes*n.nbins)
kappabin = np.zeros(n.nmodes*n.nbins) # just externally mixed aerosol,  
molwbin = np.zeros(n.nmodes*n.nbins)  # need to add mass fraction term for 
                                      # for internally mixed aerosol
for i in range(n.nmodes):
    rhobin[i*n.nbins:n.nbins+n.nbins*i] = n.rhoa[i]
    kappabin[i*n.nbins:n.nbins+n.nbins*i] = n.k[i]
    molwbin[i*n.nbins:n.nbins+n.nbins*i] = n.molw_aer[i]
    
############################### set-up grid ###################################
GRID = s.setup_grid(rhobin, kappabin, n.rkm, n.Dlow, n.nbins, 
                    n.nmodes, n.RH, n.sig, n.NAER, n.D_AER, n.T)

Y[n.nbins*n.nmodes:-3] = GRID[0] # number of aerosol in each bin
Y[:IND1] = GRID[1]               # water mass in each bin
Y_AER[0,:] = GRID[2]                  # mass of aerosol in each bin
YICE_AER[0,:] = GRID[2]                 # mass of aerosol in ice bins 

###############################################################################
def kk02(NW):
    """ Kappa Koehler theory, Petters and Kriedenwies (2007)
        
    """
    MWAT = NW*c.mw
    RH_ACT = 0
    mass_bin_centre = brac
    T = output[idx,ITEMP]
    mult = -1.0
    
    RHOAT = MWAT/c.rhow+(brac/rhobin3)
    RHOAT = (MWAT+(brac))/RHOAT

    Dw = ((MWAT + (mass_bin_centre))*6/(np.pi*RHOAT))**(1/3)
    Dd = ((mass_bin_centre*6)/(rhobin3*np.pi))**(1/3)
    KAPPA = (mass_bin_centre/rhobin3*kappabin3)/(mass_bin_centre/rhobin3)
    
    sigma = f.surface_tension(T)
    RH_EQ = mult*((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                 np.exp((4*sigma*c.mw)/(c.R*T*c.rhow*Dw)))-RH_ACT
  
    return RH_EQ


def run_sim(Y,time,Y_AER1):
        
    def dy_dt_func(t,Y):
        dy_dt = np.zeros(len(Y))
       
        svp1 = f.svp_liq(Y[ITEMP])

        # saturation ratio
        SL = svp1*Y[IRH]/(Y[IPRESS]-svp1)
        SL = (SL*Y[IPRESS]/(1 + SL))/svp1

        # water vapour mixing ratio
        WV = c.eps*Y[IRH]*svp1/(Y[IPRESS] - svp1)
        
        # CHAMBER MODEL - pressure change
        dy_dt[IPRESS] = -100*n.PRESS1*n.PRESS2*np.exp(-n.PRESS2*(time+t))

# ----------------------------change in vapour content: -----------------------
        # 1. equilibruim size of particles
        KK01 = f.kk01(Y[0:n.nbins*n.nmodes], Y[ITEMP], Y_AER1, 
                      rhobin, kappabin, molwbin)
     
        Dw = KK01[2]    # wet diameter
        RHOAT = KK01[1] # density of particles inc water and aerosol mass
        RH_EQ = KK01[0] # equilibrium diameter
    
        # 2. growth rate of particles, Jacobson p455
        # rate of change of radius
        growth_rate = f.DROPGROWTHRATE(Y[ITEMP],Y[IPRESS],SL,RH_EQ,RHOAT,Dw)
        growth_rate[np.isnan(growth_rate)] = 0# get rid of nans
        growth_rate = np.where(Y[IND1:IND2] < 1e-9, 0.0, growth_rate)
       
        # 3. Mass of water condensing
        # change in mass of water per particle
        dy_dt[:n.nbins*n.nmodes] = (np.pi*RHOAT*Dw**2)* growth_rate 
        
        # 4. Change in vapour content
        # change in water vapour mixing ratio
        dwv_dt = -1*sum(Y[IND1:IND2]*dy_dt[:IND1])       
# -----------------------------------------------------------------------------
        
        # CHAMBER MODEL - change in temperature
        dy_dt[ITEMP] = -n.Temp1*n.Temp2*np.exp(-n.Temp2*(time+t))

# --------------------------------RH change------------------------------------
        dy_dt[IRH] = svp1*dwv_dt*(Y[IPRESS]-svp1)
        dy_dt[IRH] = dy_dt[IRH] + svp1*WV*dy_dt[IPRESS]
        dy_dt[IRH] = (dy_dt[IRH] - 
                     WV*Y[IPRESS]*derivative(f.svp_liq,Y[ITEMP],dx=1.0)
                     *dy_dt[ITEMP])
        dy_dt[IRH] = dy_dt[IRH]/(c.eps*svp1**2)
# -----------------------------------------------------------------------------
        
        return dy_dt
    
#--------------------- SET-UP solver ------------------------------------------
    
    y0 = Y
    t0 = 0.0
    
    #define assimulo problem
    exp_mod = Explicit_Problem(dy_dt_func,y0,t0)
    
    # define an explicit solver
    exp_sim = CVode(exp_mod)
    exp_sim.iter = 'Newton'
    exp_sim.discr = 'BDF'
    #set parameters
    tol_list = np.zeros_like(Y)
    tol_list[0:IND1] = 1e-30 # this is now different to ACPIM (1e-25)
    tol_list[IND1:IND2] = 10 # number
    tol_list[IPRESS] = 10
    tol_list[ITEMP] =1e-4
    tol_list[IRH] = 1e-8 # set tolerance for each dydt function
    exp_sim.atol = tol_list
    exp_sim.rtol = 1.0e-8
    exp_sim.inith = 0 # initial time step-size
    exp_sim.usejac = False
    exp_sim.maxncf = 100 # max number of convergence failures allowed by solver
    exp_sim.verbosity = 40
    exp_sim.maxh = 0.1
    exp_sim.dqrhomax = 10
    exp_sim.stablimdet = True    

    t_output,y_output = exp_sim.simulate(1)
    
    return y_output[-1,:], t_output[:]


def run_sim_ice(Y,YLIQ):
     
    def dy_dt_func(t,Y):
        
        dy_dt = np.zeros(len(Y))
        
        svp = f.svp_liq(Y[ITEMP])
        svp_ice = f.svp_ice(Y[ITEMP])
        
        # vapour mixing ratio
        WV = c.eps*Y[IRH_ICE]*svp/(Y[IPRESS_ICE] - svp)
        # liquid mixing ratio
        WL = sum(YLIQ[IND1:-3]*YLIQ[0:IND1])
        # ice mixing ratio
        WI = sum(Y[IND1:IND2]*Y[0:IND1])

        Cpm = c.CP + WV*c.CPV + WL*c.CPW + WI*c.CPI
        
        # RH with respect to ice
        RH_ICE = WV/(c.eps*svp_ice/
                     (Y[IPRESS_ICE]-svp_ice)) 
# ------------------------- growth rate of ice --------------------------
        RH_EQ = 1e0 # from ACPIM, FPARCELCOLD - MICROPHYSICS.f90
        
        CAP = f.CAPACITANCE01(Y[0:IND1],np.exp(Y[IND2:-3])) 

        growth_rate = f.ICEGROWTHRATE(Y[ITEMP_ICE],Y[IPRESS_ICE],RH_ICE,RH_EQ,
                                      Y[0:IND1],np.exp(Y[IND2:-3]),CAP)
        growth_rate[np.isnan(growth_rate)] = 0# get rid of nans
        growth_rate = np.where(Y[IND1:IND2] < 1e-6, 0.0, growth_rate)
        
        # Mass of water condensing 
        dy_dt[:IND1] = growth_rate 

#---------------------------aspect ratio---------------------------------------
        DELTA_RHO = c.eps*svp/(Y[IPRESS_ICE] - svp)
        DELTA_RHOI = c.eps*svp_ice/(Y[IPRESS_ICE] - svp_ice)
        DELTA_RHO = Y[IRH_ICE]*DELTA_RHO-DELTA_RHOI
        DELTA_RHO = DELTA_RHO*Y[IPRESS_ICE]/Y[ITEMP_ICE]/c.RA
        
        RHO_DEP = f.DEP_DENSITY(DELTA_RHO,Y[ITEMP_ICE])
        
        # this is the rate of change of LOG of the aspect ratio
        dy_dt[IND2:-3] = (dy_dt[0:IND1]*(
                                       (f.INHERENTGROWTH(Y[ITEMP_ICE])-1)/
                                       (f.INHERENTGROWTH(Y[ITEMP_ICE])+2))/
                                       (Y[0:IND1]*c.rhoi*RHO_DEP))
#------------------------------------------------------------------------------        
        # Change in vapour content
        dwv_dt = -1*sum(Y[IND1:IND2]*dy_dt[0:IND1])

        # change in water vapour mixing ratio
        DRI = -1*dwv_dt
        
        dy_dt[ITEMP_ICE] = 0.0#+c.LS/Cpm*DRI
        
#---------------------------RH change------------------------------------------
        
        dy_dt[IRH_ICE] = (Y[IPRESS_ICE]-svp)*svp*dwv_dt
        
        dy_dt[IRH_ICE] = (dy_dt[IRH_ICE] - 
                          WV*Y[IPRESS_ICE]*derivative(f.svp_liq,Y[ITEMP_ICE],dx=1.0)*
                          dy_dt[ITEMP_ICE])
        dy_dt[IRH_ICE] = dy_dt[IRH_ICE]/(c.eps*svp**2)
#------------------------------------------------------------------------------        
        return dy_dt
    
    #--------------------- SET-UP solver --------------------------------------
    
    y0 = Y
    t0 = 0.0
      
    #define assimulo problem
    exp_mod = Explicit_Problem(dy_dt_func,y0,t0)
    
    # define an explicit solver
    exp_sim = CVode(exp_mod)
    exp_sim.iter = 'Newton'
    exp_sim.discr = 'BDF'
    # set tolerance for each dydt function
    tol_list = np.zeros(IND1*3+3)
    tol_list[0:IND1] = 1e-25 # mass
    tol_list[IND1:IND2] = 10 # number
    tol_list[IND2:-3] = 1e-3 # aspect ratio
    tol_list[IPRESS_ICE] = 10
    tol_list[ITEMP_ICE] =1e-4
    tol_list[IRH_ICE] = 1e-8 
    exp_sim.atol = tol_list
    exp_sim.rtol = 1.0e-8
    exp_sim.inith = 1.0e-2 # initial time step-size
    exp_sim.usejac = False
    exp_sim.maxncf = 100 # max number of convergence failures allowed by solver
    exp_sim.verbosity = 40
    exp_sim.dqrhomax = 10
    exp_sim.maxh = 0.1
    exp_sim.stablimdet = True

    t_output,y_output = exp_sim.simulate(1)
    
    return y_output[-1,:]



output = np.zeros([int(t_final/n.dt),len(Y)])
dummy = np.zeros_like(Y)
ice_aer = np.zeros([int(t_final/n.dt),IND1])
total_water_mass = np.zeros(len(output))
output_ice = np.zeros([int(t_final/n.dt),IND1*3+3])
output_ice2 = np.zeros_like(output_ice)
output_ice2[:,IND2:-3] = 1.0 # set aspect ratio to 1
dummy2 = np.zeros_like(Y)
dummy4 = np.zeros_like(Y)
ACT_DROPS = np.zeros([int(t_final/n.dt),IND1])
test = np.zeros([int(t_final/n.dt),IND1])
act_mass = np.zeros(IND1)

for idx in range(len(output)):
    if idx == 0: 
        output[idx,:] = Y
        if output[idx,ITEMP] < 273.15:
            for i, dummy3 in enumerate(rhobin):
                rhobin3 = rhobin[i]
                kappabin3 = kappabin[i]
                brac = Y_AER[idx,i]
            
                act_mass[i] = minimize_scalar(kk02,bracket=(brac*0.1, brac*50), method='brent', tol=1e-15).x*c.mw

  
            icenuc = f.icenucleation(
                    output[idx,0:IND1],Y_AER[idx,:], output[idx,IND1:IND2], output[idx,ITEMP], 
                    output[idx,IPRESS],n.nbins, n.nmodes, rhobin, 
                    kappabin, n.ncomps, n.dt,YICE_AER[idx,:], output_ice[idx,:], IND1, IND2, act_mass, output[idx,IRH])
            
            output_ice[idx,:IND1]     = icenuc[0] # ice mass
            YICE_AER[idx,:]           = icenuc[1] # aerosol mass in ice 
            output_ice[idx,IND1:IND2] = icenuc[2] # ice number
            output[idx,IND1:IND2]     = icenuc[3] # liquid number
            
            # initialise run_sim_ice with T,P,RH from run_sim (liquid)
            output_ice[idx,ITEMP_ICE]  = output[idx,ITEMP]
            output_ice[idx,IPRESS_ICE] = output[idx,IPRESS]
            output_ice[idx,IRH_ICE]    = output[idx,IRH]
            
            # check things dont go negative
            output_ice[idx,:] = np.where(output_ice[idx,:] < 0, 
                                          1e-22,
                                          output_ice[idx,:])
            
            # grow ice ------------------------------------------------------------
            dummy2[:IND2] = output[idx,:IND2] # liquid mass and drop/aerosol number
            dummy2[-3:] = output[idx,-3:] # pressure, temperature and RH
            
            output_ice[idx,:] = run_sim_ice(output_ice[idx,:],dummy2)
            # ----------------------------------------------------------------------
            
            # check things dont go negative
            output_ice[idx,:] = np.where(output_ice[idx,:] < 0,
                                         0.0,
                                         output_ice[idx,:])
            
            # change RH and Temp due to ice formation and growth
            output[idx,IRH]   = output_ice[idx,IRH_ICE]
            output[idx,ITEMP] = output_ice[idx,ITEMP_ICE]
    else:
        # solve warm cloud ODEs
        output[idx,:], test_time = run_sim(output[idx-1,:],(idx-1), Y_AER[idx-1,:])
       
        Y_AER[idx,:] = Y_AER[idx-1,:]
# -------------------- calculate the number of activated drops ----------------------
        # 1. calculate mass for activation
        act_mass = np.zeros(IND1)
       
        for i, dummy3 in enumerate(rhobin):
            rhobin3 = rhobin[i]
            kappabin3 = kappabin[i]
            brac = Y_AER[idx,i]
            
            act_mass[i] = minimize_scalar(kk02,bracket=(brac*0.1, brac*50), method='brent', tol=1e-15).x*c.mw

        # 2. find the number of aerosol with water mass greater than mass for activation 
        ACT_DROPS[idx,:] = np.where(output[idx,:IND1] > act_mass, 
                                    output[idx,IND1:IND2], 
                                    0.0)
# ---------------------------------------------------------------------------------------      
        # calculate ice if below zero degrees
        if output[idx,ITEMP] < 273.15:
            # 1. find the number of particles that freeze
            icenuc = f.icenucleation(
                    output[idx,0:IND1],Y_AER[idx-1,:], output[idx,IND1:IND2], output[idx,ITEMP], 
                    output[idx,IPRESS],n.nbins, n.nmodes, rhobin, 
                    kappabin, n.ncomps, n.dt,YICE_AER[idx-1,:], output_ice[idx-1,:], IND1, IND2, act_mass, output[idx,IRH])
            
            output_ice[idx,:IND1]     = icenuc[0] # ice mass
            YICE_AER[idx,:]           = icenuc[1] # aerosol mass in ice 
            output_ice[idx,IND1:IND2] = icenuc[2] # ice number
            output[idx,IND1:IND2]     = icenuc[3] # liquid number
            Y_AER[idx,:] = Y_AER[idx-1,:]
            
            
            # initialise run_sim_ice with T,P,RH from run_sim (liquid)
            output_ice[idx,ITEMP_ICE]  = output[idx,ITEMP]
            output_ice[idx,IPRESS_ICE] = output[idx,IPRESS]
            output_ice[idx,IRH_ICE]    = output[idx,IRH]
            
            # check things dont go negative
            output_ice[idx,:] = np.where(output_ice[idx,:] < 0, 
                                          1e-22,
                                          output_ice[idx,:])
     
            dummy2[:IND2] = output[idx,:IND2] # liquid mass and drop/aerosol number
            dummy2[-3:] = output[idx,-3:] # pressure, temperature and RH
            
            # 2. grow ice by solving cold ODEs
            output_ice[idx,:] = run_sim_ice(output_ice[idx,:],dummy2)
            
            # check things dont go negative
            output_ice[idx,:] = np.where(output_ice[idx,:] < 0,
                                         0.0,
                                         output_ice[idx,:])
            
            # change RH and Temp due to ice formation and growth
            output[idx,IRH]   = output_ice[idx,IRH_ICE]
            output[idx,ITEMP] = output_ice[idx,ITEMP_ICE]
            
            # melting -> move melted water to liquid - TO DO
            
            print('time is = ',idx)
            
        # check conserve mass, total_water_mass should be constant   
        total_water_mass[idx] = (sum(output_ice[idx,0:IND1]*output_ice[idx,IND1:IND2])+
                                 sum(output[idx,0:IND1]*output[idx,IND1:IND2])+(
                                 c.eps*output[idx,IRH]*f.svp_liq(output[idx,ITEMP])/(output[idx,IPRESS]
                                 - f.svp_liq(output[idx,ITEMP])*output[idx,IRH]))/1000)# think vapour was in g/kg, /1000 now in kg/kg, but what is mass of bins units? - 04/01/18
fig = plt.figure()
    
ax1 = fig.add_subplot(221)
ax1.plot(np.sum(ACT_DROPS, axis=1)/1e6,label='py')
#ax1.plot(nc['TIMES'],nc['NDROP'][:,0,0,0]/1e6, label = 'ACPIM')
ax1.set_title('Activated Drops')
ax1.set_ylabel('# cm$^{-3}$')
ax1.legend()
# 
ax2 = fig.add_subplot(222)
ax2.plot(np.sum(output_ice[:,IND1:IND2], axis=1)/1e6)
#ax2.plot(nc['TIMES'],nc['NICE'][:,0,0,0]/1e6)
# ax2.plot(obs)
ax2.set_title('Ice Number Concentration')
ax1.set_ylabel('# cm$^{-3}$')
# 
ax3 = fig.add_subplot(223)
ax3.plot(output[:,IRH])
#ax3.plot(nc['TIMES'],nc['RH'][:,0,0,0])
ax3.set_title('RH')
ax3.set_xlabel('Time (s)')
# 
ax4 = fig.add_subplot(224)
ax4.plot(output[:,ITEMP]-273.15)
#ax4.plot(nc['TIMES'],nc['TEMP'][:,0,0,0]-273.15)
ax4.set_title('Temperature')
ax4.set_xlabel('Time (s)')
# 
plt.savefig('output.pdf')
plt.tight_layout()
plt.show()
# =============================================================================


# =============================================================================
############################### NOTES #########################################
# 29/09/2017 ------------------------------------------------------------------
# BUG - model not running for T -8 -> 1C and RH = 0.9 in NAMELIST
# does not get past first time step
# higher RH does run
# setting growthrate to 1e-11 makes model run
# calculation of growthrate looks fine as does dy_dt[RH]
# Solver passing in NANs for temperature at some points?
# is it because growthrate goes negative?
# 25/10/2017 ------------------------------------------------------------------
# FIX - calculation of aspect ratio (instead of just setting to 1)
# model now runs for T -7 RH = 0.9 
# BUG - losing water mass somewhere when there is ice  
# BUG - Mass required for activation too high - 10/11/2017 ------------------
# 24/11/2017 ------------------------------------------------------------------
# FIX - Mass required for activation corrected, things activate, similar to
#       ACPIM results. Changed RHOAT in KK02 to RHOW and bracketed root in brent
# BUG - Is RHOAT calculated correctly? It shouldnt be too different than RHOW
# 27/11/2017 ------------------------------------------------------------------
# BUG - Freezing around an order of magnitude more ice than ACPIM?
# 14/02/2018 ------------------------------------------------------------------
# FIX - Ice number now agrees with ACPIM, changes to f.icenucleation and 
#        f.activesites made.
# BUG - won't run for small diameters and low kappa values eg 70nm and 0.0061
# FIX - changed tolerence on mass from 1e-25 to 1e-30
# 06/04/2018 ------------------------------------------------------------------
# BUG - too much ice and too high RH at low temperatures (< -25C)
# FIX - moved calculation of capacitance to inside dy_dt_func
