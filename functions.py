#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 15:04:49 2017

@author: mbexkes3
"""
import numpy as np
import constants as c
from math import sqrt

#import namelist as n

def svp_liq(T):
    """satuation vapour pressure, Buck (1996)"""
    return 100*6.1121*np.exp(((18.678-(T-273)/234.5)*((T-273)/(257.14+(T-273)))))

def svp_ice(T):
    """saturation vapour pressure over ice """
    return 100*6.1115e0*np.exp((23.036e0
        - (T-273.15)/333.7e0)*(T-273.15)/(279.82e0 
          + (T-273.15)))

def Diff(T,P):
    """Diffusion of water vapour in air"""
    T1 = max(T,200)
    return 2.11e-5*(T1/273.15)**1.94*(101325/P)
    
def KA(T):
    """Thermal conductivity of air"""
    T1 = max(T,200)
    return (5.69+0.017*(T1-273.15))*1e-3*c.JOULES_IN_A_CAL

def surface_tension(T):
    """surface tension of water - pruppacher and klett p130"""
    TC = T - 273.15
    TC = max(TC, -40)
    surface_tension = (75.93 + 0.115*TC + 6.818e-2*TC**2 +
                       6.511e-3*TC**3 + 2.933e-4*TC**4 +
                       6.283e-6*TC**5 + 5.285e-8*TC**6)
    if TC > 0:
        surface_tension = 76.1 - 0.155*TC
        
    surface_tension = surface_tension*c.JOULES_IN_AN_ERG # convert to J/cm2
    surface_tension = surface_tension*1e4 # convert to J/m2
    
    return surface_tension

def DROPGROWTHRATE(T,P,RH,RH_EQ,RHOAT,D):
    """ Jacobson 
    """
    RAD = D/2.
    RHOA = P/c.RA/T
    D1 = Diff(T,P)
    K1 = KA(T)
    FV = 1
    FH = 1
    
    SVP = svp_liq(T) # so only calls function one time
    #############
    #place to add ventilation stuff(MICROPHYSICS line 817)
    #############
    DSTAR = D1*FV/(RAD/(RAD+0.7*8e-8)+D1*FV/RAD/c.ALPHA_COND*sqrt(2*np.pi/c.RV/T))
    KSTAR = K1*FH/(RAD/(RAD+2.16e-7)+K1*FH/RAD/c.ALPHA_THERM/c.CP/RHOA*sqrt(2*np.pi/c.RA/T))
    DROPGROWTHRATE = DSTAR*c.LV*RH_EQ*SVP*c.rhow/KSTAR/T*(c.LV*c.mw/T/c.R-1)
    DROPGROWTHRATE = DROPGROWTHRATE+c.rhow*c.R*T/c.mw
    
    return DSTAR*(RH-RH_EQ)*SVP/RAD/DROPGROWTHRATE

def kk01(MWAT, T, mass_bin_centre, rhobin, kappabin, molwbin):
    """ Kappa Koehler theory, Petters and Kriedenwies (2007)
        
    """
    RHOAT = MWAT/c.rhow+(mass_bin_centre/rhobin)
    RHOAT = (MWAT+(mass_bin_centre))/RHOAT
    Dw = ((MWAT + (mass_bin_centre))*6/(np.pi*RHOAT))**(1/3)
    Dd = ((mass_bin_centre*6)/(rhobin*np.pi))**(1/3)
    KAPPA = (mass_bin_centre/rhobin*kappabin)/(mass_bin_centre/rhobin)
    sigma = surface_tension(T)
    RH_EQ = ((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                 np.exp((4*sigma*c.mw)/(c.R*T*RHOAT*Dw)))

    return RH_EQ, RHOAT, Dw, Dd

def movingcentre(mwat, MWAT_CENTRE, mwat_edges, naer, nmodes, nbins):
    
    AVEMASS = np.zeros(nmodes*nbins)
    Ynew = np.zeros(nmodes*nbins*2)
    
    # calculate the average particle mass in each bin with 1 or more particles in
    # > 1 is used instead of > 0 to avoid bins with very large particle masses
    dummy = np.zeros_like(AVEMASS)
    for i, y in enumerate(zip(mwat,naer)):
        if y[1] > 1: # y[1] is naer
            dummy[i] = y[0]/y[1] # y[0] is mwat
            
    AVEMASS[:] = np.where(naer > 1.0, dummy, MWAT_CENTRE)
    
    #set-up variables
    mwat1 = np.reshape(mwat,[nmodes,nbins])
    naer1 = np.reshape(naer,[nmodes,nbins])
    AVEMASS1 = np.reshape(AVEMASS, [nmodes, nbins])
    mwat_new = np.zeros_like(mwat1)
    naer_new = np.zeros_like(mwat1)

    for j in range(nmodes):
        for i in range(nbins):
            
            idx = (np.abs(mwat_edges - AVEMASS1[j,i])).argmin()
           
            if idx > nbins-1 or idx < 0: continue
        
            mwat_new[j,idx] = mwat_new[j,idx]+mwat1[j,i]
            naer_new[j,idx] = naer_new[j,idx]+naer1[j,i]
            
    Ynew[0:nmodes*nbins] = np.reshape(mwat_new, nmodes*nbins)
    Ynew[nmodes*nbins:] = np.reshape(naer_new, nmodes*nbins) 
    
    #check things dont go negative
  #  Ynew[:] = np.where(Ynew < 1e-40, 1e-40, Ynew)      
    return Ynew
 
def KOOPNUCLEATIONRATE(AW, T, P, nbins, nmodes):
    """ function to calculate the homogeneous freezing rate following
        Koop et al (2000) 
    """
    PG = P/1e9
    K_WATER_AMB = 1.6e0
    DK_WATER_DP = -8.8e0
    K_ICE_AMB = 0.22e0
    DK_ICE_DP = -0.17e0
    
    # see Table 1 in Koop et al (2000) for these equations
    # INTEGRAL3 = eq.4 + eq.3 - eq.5 
    INTEGRAL3=( (-230.76e0 - 0.1478e0 * T + 4099.2e0 * T**(-1) +
                 48.8341e0 * np.log(T) ) * 
               (PG - 0.5e0 * (K_WATER_AMB + DK_WATER_DP * PG) * PG**2
               - (1e0/6e0) * DK_WATER_DP * PG**3) 
               - (19.43e0 - 2.2e-3 * T + 1.08e-5 * T**2 ) *
               (PG - 0.5e0 * (K_ICE_AMB + DK_ICE_DP * PG) * PG**2 -
                (1e0/6e0) * DK_ICE_DP * PG**3 ) )
               
    MUDIFF0 = 210368e0 + 131.438e0*T - 3.32373e6*T**(-1) - 41729.1e0 * np.log(T) # eq. 2
    
    DELTAAW = AW * np.exp(INTEGRAL3 / c.R / T) - np.exp(MUDIFF0 / c.R / T) # eq.
    
    LOGJ = -906.7e0 + 8502e0*DELTAAW - 26924e0*DELTAAW**2 + 29180e0*DELTAAW**3 # eq. 7
    
    return (10e0**LOGJ)*1e6 # nucleation rate in m^-3 s^-1   

def ACTIVESITES(T, MBIN2, rhobin, nbins, nmodes, ncomps):
    """ calculate number of activesites per aerosol particle
        currently only for Feldspar """
    NS = 10e0**(-0.1963e0*T + 60.2118e0)
    return (np.pi * (6e0 * MBIN2/(np.pi*rhobin))**(2/3))*NS
   
def icenucleation(MWAT, MBIN2, n_aer_bin, T, P, nbins, nmodes, rhobin, kappabin, ncomps, dt):
    """ calculate number of ice crystals nucleated this time step
        homogeneous aw -> Koop et al (2000) 
        heterogeneous ns -> Connolly et al (2012) """
        
    YICE = np.zeros(nbins*nmodes)
    
    # calculate water activity
    RHOAT = MWAT/c.rhow+(MBIN2/rhobin)
    RHOAT = (MWAT+(MBIN2))/RHOAT
    Dw = ((MWAT + (MBIN2))*6/(np.pi*RHOAT))**(1/3)
    Dd = ((MBIN2*6)/(rhobin*np.pi))**(1/3)
    KAPPA = (MBIN2/rhobin*kappabin)/(MBIN2/rhobin)
    AW = (Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))
    
    # calculate homogeneous freezing rate following Koop et al (2000)
    # JW units m^-3 s^-1
    JW = KOOPNUCLEATIONRATE(AW, T, P, nbins, nmodes)
    
    # find the number of homogeneously frozen drops 
    #P=1-exp(-J*V*t), right hand column page 3 Koop et al
    number_frozen = n_aer_bin * (1 - np.exp(-JW * (MWAT/c.rhow) * dt))
    
    # calculate number of active sites - needs to be for each composition, just now for feldspar only
    NS = ACTIVESITES(T, MBIN2, rhobin, nbins, nmodes, ncomps)
 
    # HETEROGENEOUS CRITERIA
    NS = np.where(MWAT < 1e-22, 0.0, NS)
    
    # number of ice in each bin
    YICE[:] = number_frozen + (n_aer_bin-number_frozen)*(1e0-np.exp(-NS))
     
    return YICE

def CAPACITANCE01(MWAT, PHI):
    RHOICE = c.rhoi
    VOL = np.where(MWAT > 0, MWAT/RHOICE, 1e-50)
    
    A = (3e0*VOL/(4e0*np.pi*PHI))**(1e0/3e0)
    C = A*PHI
  
#### IF slower than everything done in a WHERE ###############################
  #  CAP = np.zeros_like(PHI)
    
   # for i in PHI:
    #    if i < 1.0:
     #       CAP = A*np.sqrt(1-i**2)/np.arcsin(np.sqrt(1-i**2))
      #  else:
       #     CAP =  C*np.sqrt(1-i**-2)/np.log((1+np.sqrt(1-i**-2))*i)
###############################################################################
       
    CAP = (np.where(PHI < 1.0, 
                    A*np.sqrt(1-PHI**2)/np.arcsin(np.sqrt(1-PHI**2)),
                    C*np.sqrt(1-PHI**-2)/np.log((1+np.sqrt(1-PHI**-2))*PHI))
                    )
    
    CAP = np.where(np.abs(PHI-1) < 1e-4,A,CAP)
    
    return CAP

def ICEGROWTHRATE(T, P, RH_ICE, RH_EQ, MWAT, AR, RAD):
    """ Jacobson, Second edition (2005) 'Fundementals of Atmos. Modelling' 
        eq.16.76, p549  """
    
    RHOICE = 910.0e0# this is not calculated yet, this is density of ice at 0C
 #   RAD = CAPACITANCE01(MWAT, RHOICE, AR)
    
    RHOA = P/c.RA/T # density of air
    D1 = Diff(T,P) # diffusivity of water vaour in air
    K1 = KA(T) # thermal conductivity of air
    FV = 1e0
    FH = 1e0
    
    #############
    #place to add ventilation stuff(MICROPHYSICS line 817)
    #############
    
    DSTAR = D1*FV/(RAD/(RAD+0.7*8e-8)+D1*FV/RAD/c.ALPHA_DEP*np.sqrt(2*np.pi/c.RV/T))
    KSTAR = K1*FH/(RAD/(RAD+2.16e-7)+K1*FH/RAD/c.ALPHA_THERM_ICE/c.CP/RHOA*np.sqrt(2*np.pi/c.RA/T))
    
    ICEGROWTHRATE = DSTAR*c.LS*RH_EQ*svp_ice(T)/KSTAR/T*(c.LS*c.mw/T/c.R-1)
    
    ICEGROWTHRATE = ICEGROWTHRATE+c.R*T/c.mw
    
    return 4e0*np.pi*RAD*DSTAR*(RH_ICE-RH_EQ)*svp_ice(T)/ICEGROWTHRATE

def INHERENTGROWTH(T):
    """ look up table from fig 12 Chen and Lamb 1994
        see CONSTANTS_VARS.f90 ACPIM GAMMA_XTAL
    """
    GAMMA_XTAL = [1.4, 0.7, 0.55, 0.55, 0.45, 0.6, 1.4, 1.0]
    T_XTAL = [233.15, 243.15, 248.15, 253.15, 258.15, 263.15, 268.15, 273.15 ]
    idx = (np.abs(np.array(T_XTAL) - T)).argmin()
    return GAMMA_XTAL[idx]

def DEP_DENSITY(DEL_RHO, T):
    DEP_D = 910.0*np.exp(-3*max(1e3*DEL_RHO - 0.05, 0)/INHERENTGROWTH(T))
    DEP_D = max(DEP_D,50)
    return DEP_D

def aspect_ratio(T, P, RH, MICE, MICEOLD, AR, QV):
    """ calculate aspect ratio, this is PHI in capacitance function,
        see MICROPHYSIC.f90 subroutine icenucleation
        and lines 1888 - 1912 """
        
    DELTA_RHO = c.eps*svp_liq(T)/(P - svp_liq(T))
    DELTA_RHOI = c.eps*svp_ice(T)/(P - svp_ice(T))
    DELTA_RHO = RH*DELTA_RHO-DELTA_RHOI
    DELTA_RHO = DELTA_RHO*P/T/c.RA
    
    RHO_DEP = DEP_DENSITY(DELTA_RHO,T)
    print(RHO_DEP)
    GAMMAICE = INHERENTGROWTH(T) #this is working
    return  AR * np.exp((GAMMAICE - 1)/(GAMMAICE + 2)*np.log(((MICE - MICEOLD)/
                                       RHO_DEP)))
   # AR = AR * np.exp((GAMMAICE - 1)/(GAMMAICE + 2)*np.log((QV + (MICE - MICEOLD)/
    #                                   RHO_DEP)/QV))
    
    
    
   
   
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    