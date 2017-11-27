#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 14:51:12 2017

@author: mbexkes3
"""
import numpy as np
from scipy.misc import derivative
from scipy.optimize import minimize
from scipy.optimize import brent

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

############# SET-UP INDEXS ###################################################
IND1 = n.nbins*n.nmodes
IND2 = IND1*2
IPRESS = -3#0 + n.nbins*n.nmodes*2
ITEMP = -2#1 + n.nbins*n.nmodes*2
IRH = -1#2 + n.nbins*n.nmodes*2
IAR = 0
IQV = 1

IPRESS_ICE = -3
ITEMP_ICE = -2
IRH_ICE = -1
###############################################################################

Y = np.zeros(n.nbins*n.nmodes*2+3)

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
MWAT_CENTRE =GRID[2] # this changes after rebinning
n_aer_bin = GRID[0] # the changes after rebinning
AER_MASS_CENTRE = GRID[5] # does not change as no semi-vols or coagulation/collisions
AER_MASS_EDGES = GRID[4] # does not change as no semi-vols or coagulation/collisions

Y[0:n.nbins*n.nmodes] = MWAT_CENTRE
Y[n.nbins*n.nmodes:-3] = n_aer_bin

###############################################################################
def kk02(MWAT):
    """ Kappa Koehler theory, Petters and Kriedenwies (2007)
        
    """
    mass_bin_centre = AER_MASS_CENTRE2
    T = output[idx,ITEMP]
    mult = -1.0
    
    RHOAT = MWAT/c.rhow+(mass_bin_centre/rhobin3)
    RHOAT = (MWAT+(mass_bin_centre))/RHOAT
    RHOAT = c.rhow
    Dw = ((MWAT + (mass_bin_centre))*6/(np.pi*RHOAT))**(1/3)
    Dd = ((mass_bin_centre*6)/(rhobin3*np.pi))**(1/3)
    KAPPA = (mass_bin_centre/rhobin3*kappabin3)/(mass_bin_centre/rhobin3)
    sigma = f.surface_tension(T)
    RH_EQ = mult*((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                 np.exp((4*sigma*c.mw)/(c.R*T*RHOAT*Dw)))

    return RH_EQ


def run_sim(Y,time):
        
    def dy_dt_func(t,Y):
        
        dy_dt = np.zeros(len(Y))
        
# ------temperature change due to expansion and condensation-------------------
        WV = c.eps*Y[IRH]*f.svp_liq(Y[ITEMP])/(Y[IPRESS] - f.svp_liq(Y[ITEMP]))
        WL = sum(Y[n.nbins*n.nmodes:-3]*Y[0:n.nbins*n.nmodes])
        Rm = c.RA+WV*c.RV
        Cpm = c.CP + WV*c.CPV + WL*c.CPW # just for vapour and liquid
        
        # ADIABATIC PARCEL
        # dy_dt[ITEMP] = (Rm*Y[ITEMP]*dy_dt[IPRESS])/(Cpm*Y[IPRESS])- c.LV*WV/Cpm
        # dy_dt[IPRESS] = -(Y[IPRESS]*c.g*n.w)/(Rm*Y[ITEMP])
        
        # CHAMBER MODEL
        dy_dt[ITEMP] = -n.Temp1*n.Temp2*np.exp(-n.Temp2*time)
        dy_dt[IPRESS] = -100*n.PRESS1*n.PRESS2*np.exp(-n.PRESS2*time)
# -----------------------------------------------------------------------------        
        
# ----------------------------change in vapour content: -----------------------
        # 1. equilibruim size of particles
        KK01 = f.kk01(Y[0:n.nbins*n.nmodes], Y[ITEMP], AER_MASS_CENTRE, 
                      rhobin, kappabin, molwbin)
        Dw = KK01[2]    # wet diameter
        RHOAT = KK01[1] # density of particles inc water and aerosol mass
        RH_EQ = KK01[0] # equilibruim diameter
      
        # 2. growth rate of particles, Jacobson p455
        # rate of change of radius
        growth_rate = f.DROPGROWTHRATE(Y[ITEMP],Y[IPRESS],Y[IRH],RH_EQ,RHOAT,Dw) 
        
        # 3. Mass of water condensing
        # change in mass of water per particle
        dy_dt[0:n.nbins*n.nmodes] = np.pi*RHOAT*Dw**2 * growth_rate 
       
        # 4. Change in vapour content
        # change in water vapour mixing ratio
        dwv_dt = -1*sum(Y[n.nbins*n.nmodes:-3]*dy_dt[0:n.nbins*n.nmodes]) 
# -----------------------------------------------------------------------------
        
# --------------------------------RH change------------------------------------
        dy_dt[IRH] = f.svp_liq(Y[ITEMP])*dwv_dt*(Y[IPRESS]-f.svp_liq(Y[ITEMP]))
        dy_dt[IRH] = dy_dt[IRH] + f.svp_liq(Y[ITEMP])*WV*dy_dt[IPRESS]
        dy_dt[IRH] = dy_dt[IRH] - WV*Y[IPRESS]*derivative(f.svp_liq,Y[ITEMP],dx=1.0)*dy_dt[ITEMP]
        dy_dt[IRH] = dy_dt[IRH]/(c.eps*f.svp_liq(Y[ITEMP])**2)
# -----------------------------------------------------------------------------
        
        dy_dt_func.wv = WV
        dy_dt_func.wl = WL
        dy_dt_func.rt = WV+WL
     
        return dy_dt
    
    #--------------------- SET-UP solver ------------------------------------------
    
    y0 = Y
    t0 = 0.0
    
    #define assimulo problem
    exp_mod = Explicit_Problem(dy_dt_func,y0,t0)
    
    # define an explicit solver
    exp_sim = CVode(exp_mod)
    
    #set parameters
    tol_list = np.zeros_like(Y)
    tol_list[0:n.nbins*n.nmodes] = 1e-25
    tol_list[n.nbins*n.nmodes:-3] = 10 # number
    tol_list[IPRESS] = 10
    tol_list[ITEMP] =1e-4
    tol_list[IRH] = 1e-8 # set tolerance for each dydt function
    exp_sim.atol = tol_list
    exp_sim.rtol = 1.0e-4
    exp_sim.inith = 1.0e-2 # initial time step-size
    exp_sim.usejac = False
    exp_sim.maxncf = 100 # max number of convergence failures allowed by solver
    
    t_output,y_output = exp_sim.simulate(1)
    
    return y_output[-1,:]


def run_sim_ice(Y,YLIQ,CAP):
     
    def dy_dt_func(t,Y):
        
        dy_dt = np.zeros(len(Y))
        
        WV = c.eps*Y[IRH_ICE]*f.svp_liq(Y[ITEMP_ICE])/(Y[IPRESS_ICE] - f.svp_liq(Y[ITEMP_ICE]))
        WL = sum(YLIQ[IND1:-3]*YLIQ[0:IND1])
        WI = sum(Y[IND1:IND2]*Y[0:IND1])
        Cpm = c.CP + WV*c.CPV + WL*c.CPW + WI*c.CPI
        
        RH_ICE = WV/(c.eps*f.svp_ice(Y[ITEMP_ICE])/(Y[IPRESS_ICE]-f.svp_ice(Y[ITEMP_ICE]))) #this is fine

        RH_EQ = 1e0 # from ACPIM, FPARCELCOLD - MICROPHYSICS.f90
        
        #pdb.set_trace() putting this here crashes kernel?
         
        growth_rate = f.ICEGROWTHRATE(Y[ITEMP_ICE],Y[IPRESS_ICE],RH_ICE,RH_EQ,
                                      Y[0:IND1],np.exp(Y[IND2:-3]),CAP)
  
        growth_rate = np.where(Y[n.nbins*n.nmodes:n.nbins*n.nmodes*2] < 1e-6, 0.0, growth_rate)
 
        # Mass of water condensing 
        dy_dt[0:n.nbins*n.nmodes] = growth_rate # growth_rate looks fine but solver can't solve it? time step goes too small
        
        #aspect ratio
        DELTA_RHO = c.eps*f.svp_liq(Y[ITEMP_ICE])/(Y[IPRESS_ICE] - f.svp_liq(Y[ITEMP_ICE]))
        DELTA_RHOI = c.eps*f.svp_ice(Y[ITEMP_ICE])/(Y[IPRESS_ICE] - f.svp_ice(Y[ITEMP_ICE]))
        DELTA_RHO = Y[IRH_ICE]*DELTA_RHO-DELTA_RHOI
        DELTA_RHO = DELTA_RHO*Y[IPRESS_ICE]/Y[ITEMP_ICE]/c.RA
        
        RHO_DEP = f.DEP_DENSITY(DELTA_RHO,Y[ITEMP_ICE]) # this is fine
        
        # this is the rate of change of LOG of the aspect ratio
        dy_dt[n.nbins*n.nmodes*2:-3] = dy_dt[0:n.nbins*n.nmodes]*(
                                       (f.INHERENTGROWTH(Y[ITEMP_ICE])-1)/
                                       (f.INHERENTGROWTH(Y[ITEMP_ICE])+2))/(Y[0:n.nbins*n.nmodes]*c.rhoi*RHO_DEP)
        
        # Change in vapour content
        dwv_dt = -1*sum(Y[IND1:IND2]*dy_dt[0:IND1])

        # change in water vapour mixing ratio
        DRI = -1*dwv_dt
        
        dy_dt[ITEMP_ICE] = c.LS/Cpm*DRI
        
        #RH change
        dy_dt[IRH_ICE] = (Y[IPRESS_ICE]-f.svp_liq(Y[ITEMP_ICE]))*f.svp_liq(Y[ITEMP_ICE])*dwv_dt
        
        dy_dt[IRH_ICE] = dy_dt[IRH_ICE] - WV*Y[IPRESS_ICE]*derivative(f.svp_liq,Y[ITEMP_ICE],dx=1.0)*dy_dt[ITEMP_ICE]
        dy_dt[IRH_ICE] = dy_dt[IRH_ICE]/(c.eps*f.svp_liq(Y[ITEMP_ICE])**2)
        
        dy_dt_func.wvi = WV
        dy_dt_func.wli = WL
        dy_dt_func.rti = WV+WL
     
        return dy_dt
    
    #--------------------- SET-UP solver ------------------------------------------
    
    y0 = Y
    t0 = 0.0
      
    #define assimulo problem
    exp_mod = Explicit_Problem(dy_dt_func,y0,t0)
    
    # define an explicit solver
    exp_sim = CVode(exp_mod)
    
    # set parameters
    # set tolerance for each dydt function
    tol_list = np.zeros(IND1*3+3)
    tol_list[0:IND1] = 1e-25 # mass
    tol_list[IND1:IND2] = 10 # number
    tol_list[IND2:-3] = 1e-3 # aspect ratio
    tol_list[IPRESS_ICE] = 10
    tol_list[ITEMP_ICE] =1e-4
    tol_list[IRH_ICE] = 1e-8 
    exp_sim.atol = tol_list
    exp_sim.rtol = 1.0e-3
    exp_sim.inith = 1.0e-2 # initial time step-size
    exp_sim.usejac = False
    exp_sim.maxncf = 100 # max number of convergence failures allowed by solver
    
    t_output,y_output = exp_sim.simulate(1)
    
    return y_output[-1,:]


t_final = 50
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

for idx in range(len(output)):
    if idx == 0: 
        output[idx,:] = run_sim(Y,idx)
        if output[idx,ITEMP] < 273: 
            # calculate number of ice crystals that nucleate on first time-step
            icenuc = f.icenucleation(output[idx,0:IND1],AER_MASS_CENTRE, 
                                 output[idx,IND1:IND2], output[idx,ITEMP], output[idx,IPRESS],
                                 n.nbins, n.nmodes, rhobin, kappabin, n.ncomps, n.dt)
            
            icenuc[2*n.nbins:3*n.nbins] = 0.0
            
            if output[idx,IRH] < 1.0:
                icenuc[:] = 0.0

            # ice number
            output_ice2[idx,IND1:IND2] = icenuc
            # ice mass
            output_ice2[idx,0:IND1] = output[idx,0:IND1] # output-ice mass should ne re-initailised with liquid grid, this stops ice bins growing through ice processes and makes them grow as if they were liquid bins
            
            # take new ice away from liquid, number
            output[idx,IND1:IND2] = output[idx,IND1:IND2] - icenuc 
            
            output_ice[idx,0:IND2] = output_ice2[idx,0:IND2] # number and mass
            output_ice[idx,IND2:-3] = np.log(1.0) # set initial aspect ratio to 1, need to change this 
            output_ice[idx,-3:-1] = output_ice2[idx,-3:-1] # PRESS, TEMP and RH
            
    else:
        # solve warm cloud ODEs
        output[idx,:] = run_sim(output[idx-1,:],idx)
      
        
######################## calculate the number of activated drops############################
     
        act_mass = np.zeros(IND1)
        # this needs fixing, act_mass is too high!
        for i, dummy3 in enumerate(rhobin):
            rhobin3 = rhobin[i]
            kappabin3 = kappabin[i]
            AER_MASS_CENTRE2 = AER_MASS_CENTRE[i]
            act_mass[i] = brent(kk02, tol=1e-15,brack = (AER_MASS_CENTRE2*0.5, AER_MASS_CENTRE2*50), maxiter = 100)
            
        ACT_DROPS[idx,:] = np.where(output[idx,:IND1] > act_mass, output[idx,IND1:IND2], 0.0)
        
#############################################################################################      
        
        if output[idx,ITEMP] < 273:
           
            # nulceate ice
            icenuc = f.icenucleation(output[idx,0:IND1],AER_MASS_CENTRE, 
                                 output[idx,IND1:IND2], output[idx,ITEMP], output[idx,IPRESS],
                                 n.nbins, n.nmodes, rhobin, kappabin, n.ncomps, n.dt)
          
            icenuc = np.where(output[idx,:IND1] < act_mass, 0.0, icenuc)
            icenuc[2*n.nbins:3*n.nbins] = 0.0
            
            # add newly nucleated ice onto existing ice 
            output_ice2[idx,IND1:IND2] = output_ice[idx-1,IND1:IND2]+icenuc
            output_ice2[idx,0:IND1] = output[idx,0:IND1] # output-ice mass should ne re-initailised with liquid grid, this stops ice bins growing through ice processes and makes them grow as if they were liquid bins
           
            # take new ice away from liquid
            output[idx,IND1:IND2] = output[idx,IND1:IND2] - icenuc 
            
            # initialise run_sim_ice with T,P,RH from run_sim (liquid)
            output_ice2[idx,ITEMP_ICE]  = output[idx,ITEMP]
            output_ice2[idx,IPRESS_ICE] = output[idx,IPRESS]
            output_ice2[idx,IRH_ICE]    = output[idx,IRH]
            
            # check things dont go negative
            output_ice2[idx,:] = np.where(output_ice2[idx,:] < 0, 1e-50, output_ice2[idx,:])
            
            CAP = f.CAPACITANCE01(output_ice2[idx,0:IND1],output_ice2[idx,IND2:-3])
            
            # grow ice
            dummy2[:IND2] = output[idx,:IND2]
            dummy2[-3:-1] = output[idx,-3:-1]
            output_ice[idx,:] = run_sim_ice(output_ice2[idx,:],dummy2, CAP)
            
            # check things dont go negative
            output_ice[idx,:] = np.where(output_ice[idx,:] < 0, 0.0, output_ice[idx,:])
            
            # change RH and Temp due to ice formation and growth
            output[idx,IRH]   = output_ice[idx,IRH_ICE]
            output[idx,ITEMP] = output_ice[idx,ITEMP_ICE]
            
            # melting -> move melted water to liquid
         
             
        # rebin water and ice
      #  dummy = f.movingcentre(output[idx,0:IND1], MWAT_CENTRE, 
      #                          MWATGRID, output[idx,IND1:-3], n.nmodes, n.nbins)
      #  output[idx,0:IND2] = dummy 
      #  dummy = f.movingcentre(output_ice[idx,0:IND1], MWAT_CENTRE, 
      #                         MWATGRID, output_ice[idx,IND1:IND2], n.nmodes, n.nbins)
      #  output_ice[idx,0:IND2] = dummy
        
    # check conserve mass, total_water_mass should be constant   
    total_water_mass[idx] = (sum(output_ice[idx,0:IND1]*output_ice[idx,IND1:IND2])+
                                 sum(output[idx,0:IND1]*output[idx,IND1:IND2])+
                                 c.eps*output[idx,IRH]*f.svp_liq(output[idx,ITEMP])/(output[idx,IPRESS]
                                 - f.svp_liq(output[idx,ITEMP])))
    
    
# plot output
fig = plt.figure()
    
ax1 = fig.add_subplot(221)
ax1.plot(np.sum(ACT_DROPS, axis=1)/1e6)
ax1.set_title('Activated Drops')
ax1.set_ylabel('# cm$^{-3}$')

ax2 = fig.add_subplot(222)
ax2.plot(np.sum(output_ice[:,IND1:IND2], axis=1)/1e6)
ax2.set_title('Ice Number Concentration')
ax1.set_ylabel('# cm$^{-3}$')

ax3 = fig.add_subplot(223)
ax3.plot(output[:,IRH])
ax3.set_title('RH')
ax3.set_xlabel('Time (s)')

ax4 = fig.add_subplot(224)
ax4.plot(output[:,ITEMP]-273.15)
ax4.set_title('Temperature')
ax4.set_xlabel('Time (s)')

plt.savefig('output.pdf')
plt.tight_layout()
plt.show()
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
# ISSUE - losing water mass somewhere when there is ice  
# ISSUE - Mass required for activation too high - 10/11/2017 ------------------
# 24/11/2017 ------------------------------------------------------------------
# FIX - Mass required for activation corrected, things activate, similar to
#       ACPIM results. Changed RHOAT in KK02 to RHOW and bracketed root in brent
# ISSUE - Is RHOAT calculated correctly? It shouldnt be too different than RHOW
# 27/11/2017 ------------------------------------------------------------------
# ISSUE - Freezing around an order of magnitude more ice than ACPIM?

    