#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 15:04:49 2017

@author: mbexkes3
"""
import numpy as np
import constants as c
#import namelist as n

def svp_liq(T):
    """satuation vapour pressure, Buck (1996)"""
    return 0.61121*np.exp(((18.678-(T-273)/234.5)*((T-273)/(257.14+(T-273)))))

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
    RAD = D/2.
    RHOA = P/c.RA/T
    D1 = Diff(T,P)
    K1 = KA(T)
    FV = 1
    FH = 1
    #############
    #place to add ventilation stuff(MICROPHYSICS line 817)
    #############
    DSTAR = D1*FV/(RAD/(RAD+0.7*8e-8)+D1*FV/RAD/c.ALPHA_COND*np.sqrt(2*np.pi/c.RV/T))
    KSTAR = K1*FH/(RAD/(RAD+2.16e-7)+K1*FH/RAD/c.ALPHA_THERM/c.CP/RHOA*np.sqrt(2*np.pi/c.RA/T))
    DROPGROWTHRATE = DSTAR*c.LV*RH_EQ*svp_liq(T)*RHOAT/KSTAR/T*(c.LV*c.mw/T/c.R-1)
    DROPGROWTHRATE = DROPGROWTHRATE+RHOAT*c.R*T/c.mw
    return DSTAR*(RH-RH_EQ)*svp_liq(T)/RAD/DROPGROWTHRATE

def kk01(MWAT, T, mass_bin_centre, rhobin, kappabin, molwbin):
    RHOAT = MWAT/c.rhow+(mass_bin_centre/rhobin)
    RHOAT = (MWAT+(mass_bin_centre))/RHOAT
    Dw = ((MWAT + (mass_bin_centre))*6/(np.pi*RHOAT))**(1/3)
    #Dd = ((mass_bin_centre/rhobin)*6/np.pi)**(1/3)
    Dd = ((mass_bin_centre*6)/(rhobin*np.pi))**(1/3)
    KAPPA = (mass_bin_centre/rhobin*kappabin)/(mass_bin_centre/rhobin)
    sigma = surface_tension(T)
    RH_EQ = ((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                 np.exp((4*sigma*c.mw)/(c.R*T*RHOAT*Dw)))

    return RH_EQ, RHOAT, Dw, Dd
