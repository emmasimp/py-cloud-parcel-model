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
#import pickle

import constants as c
import namelist as n
import functions as f
import setup_grids as s

import importlib # NOTE: need to be in the correct directory for this to work
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
def kk02(MWAT):
    """ Kappa Koehler theory, Petters and Kriedenwies (2007)
        
    """
    mass_bin_centre = brac
    T = output[idx,ITEMP]
    mult = -1.0
    
    RHOAT = MWAT/c.rhow+(Y_AER[idx,:]/rhobin3)
    RHOAT = (MWAT+(Y_AER[idx,:]))/RHOAT
    RHOAT = c.rhow
    Dw = ((MWAT + (mass_bin_centre))*6/(np.pi*RHOAT))**(1/3)
    Dd = ((mass_bin_centre*6)/(rhobin3*np.pi))**(1/3)
    KAPPA = (mass_bin_centre/rhobin3*kappabin3)/(mass_bin_centre/rhobin3)
    sigma = f.surface_tension(T)
    RH_EQ = mult*((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                 np.exp((4*sigma*c.mw)/(c.R*T*c.rhow*Dw)))

    return RH_EQ


def run_sim(Y,time,svp1,Y_AER):
        
    def dy_dt_func(t,Y):
        dy_dt = np.zeros(len(Y))
# ------temperature change due to expansion and condensation-------------------
        WV = c.eps*Y[IRH]*svp1/(Y[IPRESS] - svp1)
        
        # CHAMBER MODEL
        dy_dt[ITEMP] = -n.Temp1*n.Temp2*np.exp(-n.Temp2*time)
        dy_dt[IPRESS] = -100*n.PRESS1*n.PRESS2*np.exp(-n.PRESS2*time)
# -----------------------------------------------------------------------------        
        
# ----------------------------change in vapour content: -----------------------
        # 1. equilibruim size of particles
        KK01 = f.kk01(Y[0:n.nbins*n.nmodes], Y[ITEMP], Y_AER, 
                      rhobin, kappabin, molwbin)
        Dw = KK01[2]    # wet diameter
        RHOAT = KK01[1] # density of particles inc water and aerosol mass
        RH_EQ = KK01[0] # equilibruim diameter
      
        # 2. growth rate of particles, Jacobson p455
        # rate of change of radius
        growth_rate = f.DROPGROWTHRATE(Y[ITEMP],Y[IPRESS],Y[IRH],RH_EQ,RHOAT,Dw)
        growth_rate[np.isnan(growth_rate)] = 0 # get rid of nans
        growth_rate = np.where(Y[IND1:IND2] < 1e-9, 0.0, growth_rate)
     
        # 3. Mass of water condensing
        # change in mass of water per particle
        dy_dt[0:n.nbins*n.nmodes] = (np.pi*RHOAT*Dw**2)* growth_rate 
       
        # 4. Change in vapour content
        # change in water vapour mixing ratio
        dwv_dt = -1*sum(Y[n.nbins*n.nmodes:-3]*dy_dt[0:n.nbins*n.nmodes]) 
# -----------------------------------------------------------------------------
        
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
        
        svp = f.svp_liq(Y[ITEMP])
        svp_ice = f.svp_ice(Y[ITEMP])
        
        WV = c.eps*Y[IRH_ICE]*svp/(Y[IPRESS_ICE] - svp)
        WL = sum(YLIQ[IND1:-3]*YLIQ[0:IND1])
        WI = sum(Y[IND1:IND2]*Y[0:IND1])
        Cpm = c.CP + WV*c.CPV + WL*c.CPW + WI*c.CPI
        
        RH_ICE = WV/(c.eps*svp_ice/
                     (Y[IPRESS_ICE]-svp_ice)) 

        RH_EQ = 1e0 # from ACPIM, FPARCELCOLD - MICROPHYSICS.f90
                 
        growth_rate = f.ICEGROWTHRATE(Y[ITEMP_ICE],Y[IPRESS_ICE],RH_ICE,RH_EQ,
                                      Y[0:IND1],np.exp(Y[IND2:-3]),CAP)
        #growth_rate = 1e-19 # this makes ACTDROPS and RH agree with ACPIM
        growth_rate = np.where(Y[IND1:IND2] < 1e-6, 0.0, growth_rate)
        
        # Mass of water condensing 
        dy_dt[0:IND1] = growth_rate 
        
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
        
        dy_dt[ITEMP_ICE] = 0.0#c.LS/Cpm*DRI
#---------------------------RH change------------------------------------------
        dy_dt[IRH_ICE] = (Y[IPRESS_ICE]-svp)*svp*dwv_dt
        
        dy_dt[IRH_ICE] = (dy_dt[IRH_ICE] - 
                          WV*Y[IPRESS_ICE]*derivative(f.svp_liq,Y[ITEMP_ICE],dx=1.0)*
                          dy_dt[ITEMP_ICE])
        dy_dt[IRH_ICE] = dy_dt[IRH_ICE]/(c.eps*svp**2)
#------------------------------------------------------------------------------        
        dy_dt_func.wvi = WV
        dy_dt_func.wli = WL
        dy_dt_func.rti = WV+WL
        
        return dy_dt
    
    #--------------------- SET-UP solver --------------------------------------
    
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
    exp_sim.rtol = 1.0e-4
    exp_sim.inith = 1.0e-2 # initial time step-size
    exp_sim.usejac = False
    exp_sim.maxncf = 100 # max number of convergence failures allowed by solver
    
    t_output,y_output = exp_sim.simulate(1)
    
    return y_output[-1,:]


#t_final = 3
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
        svp = f.svp_liq(Y[ITEMP])
        output[idx,:] = run_sim(Y,idx,svp, Y_AER[idx,:])
        if output[idx,ITEMP] < 273.15: 
            # calculate number of ice crystals that nucleate on first time-step
            icenuc = f.icenucleation(
                    output[idx,0:IND1],Y_AER[idx,:], output[idx,IND1:IND2], output[idx,ITEMP], 
                    output[idx,IPRESS],n.nbins, n.nmodes, rhobin, 
                    kappabin, n.ncomps, n.dt,YICE_AER[idx,:], output_ice[idx,:], IND1, IND2)
            
            output_ice[idx,:IND1]     = icenuc[0] # ice mass
            YICE_AER[idx,:]           = icenuc[1] # aerosol mass in ice 
            output_ice[idx,IND1:IND2] = icenuc[2] # ice number
            output[idx,IND1:IND2]     = icenuc[3] # liquid number
            # heterogeneous freezing criteria
            ALPHA_CRIT = 1
            Dd = ((Y_AER*6)/(rhobin*np.pi))**(1/3)
            MWC = (ALPHA_CRIT/6)*np.pi*Dd**3*c.rhow
            output_ice[idx,IND1:IND2] = np.where(output[idx,:IND1] < 1e-15, 0.0, 
                                             output_ice[idx,IND1:IND2])
            
            output_ice[idx,IND2:-3] = np.log(1.0) # set initial aspect ratio to 1 
            #output_ice[idx,-3:-1] = output[idx,-3:-1]# output_ice2[idx,-3:-1] # PRESS, TEMP and RH
            output_ice[idx,IRH_ICE] = output[idx,IRH]
            output_ice[idx,IPRESS_ICE] = output[idx,IPRESS]
            output_ice[idx,ITEMP_ICE] = output[idx,ITEMP]
            
    else:
        svp = f.svp_liq(output[idx-1,ITEMP])
        # solve warm cloud ODEs
        output[idx,:] = run_sim(output[idx-1,:],idx, svp, Y_AER[idx-1,:])
        Y_AER[idx,:] = Y_AER[idx-1,:]
######################## calculate the number of activated drops ##############
        act_mass = np.zeros(IND1)
       
        for i, dummy3 in enumerate(rhobin):
            rhobin3 = rhobin[i]
            kappabin3 = kappabin[i]
            brac = Y_AER[idx,i]
            act_mass[i] = brent(kk02, tol=1e-15,brack = (brac*0.5, 
                                                         brac*50), 
                                                         maxiter = 100)
            
        ACT_DROPS[idx,:] = np.where(output[idx,:IND1] > act_mass, 
                                    output[idx,IND1:IND2], 
                                    0.0)
###############################################################################      
        
        if output[idx,ITEMP] < 273.15:
               
            icenuc = f.icenucleation(
                    output[idx,0:IND1],Y_AER[idx-1,:], output[idx,IND1:IND2], output[idx,ITEMP], 
                    output[idx,IPRESS],n.nbins, n.nmodes, rhobin, 
                    kappabin, n.ncomps, n.dt,YICE_AER[idx-1,:], output_ice[idx-1,:], IND1, IND2)
            
            output_ice[idx,:IND1]     = icenuc[0] # ice mass
            YICE_AER[idx,:]           = icenuc[1] # aerosol mass in ice 
            output_ice[idx,IND1:IND2] = icenuc[2] # ice number
            output[idx,IND1:IND2]     = icenuc[3] # liquid number
            Y_AER[idx,:] = Y_AER[idx-1,:]
            
            output_ice[idx,IRH_ICE] = output[idx,IRH]
            output_ice[idx,IPRESS_ICE] = output[idx,IPRESS]
            output_ice[idx,ITEMP_ICE] = output[idx,ITEMP]
            
            # HET CRITERIA ###################################################
            # thredhold water mass
            ALPHA_CRIT = 1
            Dd = ((Y_AER*6)/(rhobin*np.pi))**(1/3)
            MWC = (ALPHA_CRIT/6)*np.pi*Dd**3*c.rhow
            # only activated drops, use act_mass variable
          #  output_ice[idx,IND1:IND2] = np.where(output[idx,:IND1] < act_mass,
           #                             0.0,
           #                             output_ice[idx,IND1:IND2])
            
            # calculate new aspect ratio of each bin with new ice in
            output_ice[idx,IND2:-3] = np.where(
                                            output_ice[idx-1, IND1:IND2]+output_ice[idx,IND1:IND2] > 1.0,
                                            (output_ice[idx-1,IND1:IND2]*output_ice[idx-1,IND2:-3]
                                            + output_ice[idx,IND1:IND2]*1.0) / output_ice[idx,IND1:IND2],
                                            1.0)
           
              
            # initialise run_sim_ice with T,P,RH from run_sim (liquid)
            output_ice[idx,ITEMP_ICE]  = output[idx,ITEMP]
            output_ice[idx,IPRESS_ICE] = output[idx,IPRESS]
            output_ice[idx,IRH_ICE]    = output[idx,IRH]
            
            # check things dont go negative
            output_ice[idx,:] = np.where(output_ice[idx,:] < 0, 
                                          1e-50,
                                          output_ice[idx,:])
            # call capacitance function with ice mass and PHI
            CAP = f.CAPACITANCE01(output_ice[idx,0:IND1],np.exp(output_ice[idx,IND2:-3]))
            
            # grow ice
            dummy2[:IND2] = output[idx,:IND2] # ice mass and ice number
            dummy2[-3:] = output[idx,-3:] # pressure, temperature and RH
            
            output_ice[idx,:] = run_sim_ice(output_ice[idx,:],dummy2, CAP)
            
            # check things dont go negative
            output_ice[idx,:] = np.where(output_ice[idx,:] < 0,
                                         0.0,
                                         output_ice[idx,:])
            
            # change RH and Temp due to ice formation and growth
            output[idx,IRH]   = output_ice[idx,IRH_ICE]
            output[idx,ITEMP] = output_ice[idx,ITEMP_ICE]
            
            # melting -> move melted water to liquid
         
      # check conserve mass, total_water_mass should be constant   
        total_water_mass[idx] = (sum(output_ice[idx,0:IND1]*output_ice[idx,IND1:IND2])+
                                 sum(output[idx,0:IND1]*output[idx,IND1:IND2])+(
                                 c.eps*output[idx,IRH]*f.svp_liq(output[idx,ITEMP])/(output[idx,IPRESS]
                                 - f.svp_liq(output[idx,ITEMP])*output[idx,IRH]))/1000)# think vapour was in g/kg, /1000 now in kg/kg, but what is mass of bins units? - 04/01/18
# =============================================================================
# pickle.dump(output, open("output01.p", "wb"))
# pickle.dump(output_ice, open("output01_ice.p", "wb"))
# 
# pickle.dump(ACT_DROPS, open("ACT_DROPS", "wb"))
# =============================================================================
    
# plot output
# =============================================================================
# =============================================================================
# fig = plt.figure()
# #     
# ax1 = fig.add_subplot(221)
# ax1.plot(np.sum(ACT_DROPS, axis=1)/1e6)
# ax1.plot(nc['TIMES'],nc['NDROP'][:,0,0,0]/1e6)
# ax1.set_title('Activated Drops')
# ax1.set_ylabel('# cm$^{-3}$')
# # 
# ax2 = fig.add_subplot(222)
# ax2.plot(np.sum(output_ice[:,IND1:IND2], axis=1)/1e6)
# ax2.plot(nc['TIMES'],nc['NICE'][:,0,0,0]/1e6)
# # ax2.plot(obs)
# ax2.set_title('Ice Number Concentration')
# ax1.set_ylabel('# cm$^{-3}$')
# # 
# ax3 = fig.add_subplot(223)
# ax3.plot(output[:,IRH])
# ax3.plot(nc['TIMES'],nc['RH'][:,0,0,0])
# ax3.set_title('RH')
# ax3.set_xlabel('Time (s)')
# # 
# ax4 = fig.add_subplot(224)
# ax4.plot(output[:,ITEMP]-273.15)
# ax4.plot(nc['TIMES'],nc['TEMP'][:,0,0,0]-273.15)
# ax4.set_title('Temperature')
# ax4.set_xlabel('Time (s)')
# # 
# # plt.savefig('output.pdf')
# plt.tight_layout()
# plt.show()
# 
# =============================================================================
fig = plt.figure()
ax1 = fig.add_subplot(221)
ax1.plot(output[:,IRH])
ax1.set_title('RH')

ax2 = fig.add_subplot(222)
ax2.plot(total_water_mass)
ax2.set_title('total water mass')

ax3 = fig.add_subplot(223)
ax3.plot(np.sum(output[:,:IND1]*output[:,IND1:IND2], axis=1))
ax3.set_title('water mass')

ax4 = fig.add_subplot(224)
ax4.plot(np.sum(output_ice[:,:IND1]*output_ice[:,IND1:IND2], axis=1))
ax4.set_title('ice mass')
plt.tight_layout()
plt.savefig('mass_compare.pdf')
plt.show()

# =============================================================================
# fig = plt.figure()
# ax1 = fig.add_subplot(221)
# ax1.plot(nc['NBIN1'][28,:,0,0,0,0],label='28')
# ax1.plot(nc['NBIN1'][29,:,0,0,0,0], label='29')
# ax1.plot(nc['NBIN1'][30,:,0,0,0,0], label='30')
# ax1.legend()
# ax1.set_title('ACPIM liquid')
# 
# ax2 = fig.add_subplot(222)
# ax2.plot(output[25,IND1:IND1+70], label='26')
# ax2.plot(output[26,IND1:IND1+70], label='27')
# ax2.plot(output[27,IND1:IND1+70], label='28')
# ax2.legend()
# ax2.set_title('pyACPIM liquid')
# 
# ax3 = fig.add_subplot(223)
# ax3.plot(nc['NBIN2'][28,:,0,0,0,0],label='28')
# ax3.plot(nc['NBIN2'][29,:,0,0,0,0], label='29')
# ax3.plot(nc['NBIN2'][30,:,0,0,0,0], label='30')
# ax3.legend()
# ax3.set_title('ACPIM ice')
# 
# ax4 = fig.add_subplot(224)
# ax4.plot(output_ice[25,IND1:IND1+70], label='26')
# ax4.plot(output_ice[26,IND1:IND1+70], label='27')
# ax4.plot(output_ice[27,IND1:IND1+70], label='28')
# ax4.legend()
# ax4.set_title('pyACPIM ice')
# plt.tight_layout()
# plt.savefig('bin_number_compare.pdf')
# plt.show()
# =============================================================================

# 
# =============================================================================
#  fig, ax1 = plt.subplots()
# # 
#  ax2 = ax1.twinx()
#  ax1.plot(np.sum(output[:,:IND1]*output[:,IND1:IND2],axis=1), label='total water')
#  ax1.plot(np.sum(output_ice[:,IND1:IND2],axis=1)/1e6, label='Ice crystals')
# # 
#  ax2.plot(output[:,IRH], '-r',label='RH')
# # ax2.set_xlim([140, 170])
# # 
#  plt.show()
# # 
#  fig, ax1 = plt.subplots()
# # 
#  ax2 = ax1.twinx()
#  ax1.plot(total_water_mass)
#  ax2.plot(np.sum(output_ice[:,IND1:IND2],axis=1)/1e6, label='Ice crystals')
# # 
#  #ax2.plot(output[:,IRH], '-r',label='RH')
#  ax1.legend()
# # 
#  plt.show()
# # 
# =============================================================================
# plt.plot(output_ice[20, IND1+n.nbins:IND2], label='1')
# plt.plot(output_ice[25, IND1+n.nbins:IND2], label = '2')
# plt.plot(output_ice[30, IND1+n.nbins:IND2], label = '3' )
# plt.plot(output_ice[35, IND1+n.nbins:IND2], label = '4')
# plt.plot(output_ice[40, IND1+n.nbins:IND2], label='5')
# plt.plot(output_ice[45, IND1+n.nbins:IND2], label= '6')
# plt.plot(output_ice[49, IND1+n.nbins:IND2], label = '7')
# plt.legend()
# plt.show()
# 
 #fig, ax1 = plt.subplots()
 #ax1.plot(np.sum(output[:,IND1:IND2],axis=1)/1e6, 'r')
 #ax2 = ax1.twinx()
 #ax2.plot(np.sum(ACT_DROPS,axis=1)/1e6)
 #plt.show()
# 
# fig, ax1 = plt.subplots()
# ax1.plot(nc['TIMES'],nc['NICE'][:,0,0,0], 'r')
# ax2 = ax1.twinx()
# ax2.plot(nc['TIMES'],nc['NDROP'][:,0,0,0])
# plt.show()
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
    
