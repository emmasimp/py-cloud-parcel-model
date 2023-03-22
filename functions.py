#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 25 15:04:49 2017

@author: mbexkes3
"""
import numpy as np
import constants as c
from math import sqrt
import variables as v

import namelist as n 
from scipy.optimize import minimize_scalar
import bisect

def svp_liq(T):
    """satuation vapour pressure, Buck (1996)"""
    return 100*6.1121*np.exp(((18.678-(T-273.15)/234.5)*((T-273.15)/(257.14+(T-273.15)))))

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

def VISCOSITY_AIR(T):
    TC = T - 273.15
    TC = max(TC, -200.0)
    
    if TC >= 0.0:
        V = (1.7180 + 0.0049*TC)*1e-5
    else:
        V = (1.7180 + 0.0049*TC-1.2e-5 *TC**2)*1e-5
    return V    

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

def VENTILATION01(DIAM, RHOAT, T, P, nbins, nmodes):
    # density of air
    RHOA = P/c.RA/T
    # diffusivity of water vapour in air
    D1 = Diff(T, P)
    # conductivity of air
    K1 = KA(T)
    # Viscosity of air
    ETA = VISCOSITY_AIR(T)    
    # Kinematic viscosity
    NU = ETA / RHOA
    # Schmitt Numbers
    NSc1 = NU / D1
    NSc2 = NU / K1
    
    ################ terminal velocity of water drops #########################
    DIAM2 = DIAM # temp array that can be changed
    NRE = np.zeros(nbins*nmodes)
    
    DIAM2 = np.where(DIAM2 > 7000e-6,
                     7000e-6,
                     DIAM2)
    MASS = np.pi/6*DIAM2**3*RHOAT
    SIGMA = surface_tension(T)
    
    # Regime 3: equations 5-12, 10-146 and 10-148 from P+K
    PHYSNUM = (SIGMA**3)*(RHOA**2)/((ETA**4)*c.g*(c.rhow-RHOA))
    PHYS6 = PHYSNUM**(1/6)
    
    for d in DIAM2:
        if d > 1070e-6:
            BONDNUM = (4.0/3.0)*c.g*(c.rhow - RHOA)*(d**2.0)/SIGMA
            X = np.log(BONDNUM*PHYS6)
            Y = (-5.00015 + 5.23778*X - 2.04914*X*X+0.475294*(X**3)
                 - 0.542819e-1*(X**4.0) + 0.23844e-2*(X**5.0))
            NRE = PHYS6*np.exp(Y)
            VEL = ETA * NRE/(RHOA*d)
            CD = 8.0*MASS*c.g*RHOA/(np.pi*((d/2.0)*ETA)**2.0)
            CD = CD / (NRE**2.0)
        if d < 1070e-6 and d > 20e-6:
            BESTNM = 32.0*((d/2.0)**3.0)*(c.rhow-RHOA)*RHOA*c.g/(3.0*ETA**2.0)
            X = np.log(BESTNM)
            Y = (-3.18657+0.992696*X-0.153193e-2*X*X
                 -0.987059e-3*(X**3.0) - 0.578878e-3*(X**4)
                 + 0.855176e-4*(X**5) - 0.327815e-5 * (X**6.0))
            NRE = np.exp(Y)
            VEL = ETA * NRE/(2.0*RHOA*(d/2.0))
            CD = BESTNM/(NRE**2.0)
        if d < 20e-6:
            MFPATH = 6.6e-8 * (101325/P)*(T/293.15)
            VEL = 2.0*((d/2.0)**2.0)*c.g*(c.rhow-RHOA)/(9.0*ETA)
            VEL = VEL * (1.0 + 1.26 * MFPATH/ (d/2.0))
            NRE = VEL * RHOA * d/ETA
            CD = 8.0*MASS*c.g*RHOA/(np.pi*((d/2.0)*ETA)**2.0)
            CD = CD/(NRE**2.0)
            
    VEL = np.where(np.isnan(VEL), 0.0, VEL) # replace nans with zeros
  #############################################################################
    
    CALC = (NSc1**(1.0/3.0))*np.sqrt(NRE)
    CALC = np.where(CALC > 51.4, 51.4, CALC)
    FV = np.where(CALC < 1.4, 
                  1.0+0.108*CALC**2,
                  0.78+0.308*CALC)
  
    CALC = (NSc2**(1.0/3.0)) * np.sqrt(NRE)
    CALC = np.where(CALC > 51.4, 51.4, CALC)
    FH = np.where(CALC < 1.4,
                  1.00+0.108*CALC**2,
                  0.78+0.308*CALC)
    return FH, FV

#def VENTILATION02(MWAT, T, P, PHI, RHOI,nbins, nmodes ) # need to add ventilation of ice
    # crystals, requires many functions and new variables to be calculated. 
            
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
    
   # FH, FV = VENTILATION01(D, RHOAT, T, P)
      
    DSTAR = D1*FV/(RAD/(RAD+0.7*8e-8)+D1*FV/RAD/c.ALPHA_COND*sqrt(2*np.pi/c.RV/T))
    
    KSTAR = K1*FH/(RAD/(RAD+2.16e-7)+K1*FH/RAD/c.ALPHA_THERM/c.CP/RHOA*sqrt(2*np.pi/c.RA/T))
    
    DROPGROWTHRATE = DSTAR*c.LV*RH_EQ*SVP*c.rhow/KSTAR/T*(c.LV*c.mw/T/c.R-1)
    
    DROPGROWTHRATE = DROPGROWTHRATE+c.rhow*c.R*T/c.mw
    
    return DSTAR*(RH-RH_EQ)*SVP/RAD/DROPGROWTHRATE

def kk01(MWAT, T, mass_bin_centre, rhobin, kappabin):
    """ Kappa Koehler theory, Petters and Kriedenwies (2007)
        
    """
    
    RHOAT = MWAT/c.rhow+(mass_bin_centre/rhobin)
    RHOAT = (MWAT+(mass_bin_centre))/RHOAT
   # RHOAT = c.rhow
    #RHOAT = 1690 
    
    Dw = ((MWAT + (mass_bin_centre))*6/(np.pi*RHOAT))**(1/3) # nearly the same
    Dd = ((mass_bin_centre*6)/(rhobin*np.pi))**(1/3) # exactly the same
    #KAPPA = (mass_bin_centre/rhobin*kappabin)/(mass_bin_centre/rhobin)
    KAPPA = kappabin
    sigma = surface_tension(T) # 
    
    RH_EQ = ((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                 np.exp((4*sigma*c.mw)/(c.R*T*c.rhow*Dw)))

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

def ACTIVESITES(T, MBIN2, rhobin, nbins, nmodes, ncomps,dt,Mass_frac):
    """ calculate number of activesites per aerosol particle
        based on ns parameterisation, given INP type in namelist
        and dictionary of coded INP types in constants.py"""
    
    NS1 = np.zeros(ncomps)
    ACT_SITE = np.zeros([nmodes,nbins])
    NS_out = np.zeros(nbins*nmodes)
    M = np.reshape(MBIN2,[nmodes,nbins])
    
    for mode in range(nmodes):    
        for i, INP_type in enumerate(c.aerosol_dict.keys()):

             # if the mode is not an INP, NS = 0
             if INP_type not in c.ns_dict:
                 NS1[i] = 0.0
       
             # select the ns equation based on first value in ns_dict
             # Equation 5 from Niemand et al (2012)
             elif c.ns_dict[INP_type][0] == 0:
                 NS1[i] = np.exp(c.ns_dict[INP_type][1]*(T-273.15)+c.ns_dict[INP_type][2])
                 
             elif c.ns_dict[INP_type][0] == 1:
                 NS1[i] = np.exp(c.ns_dict[INP_type][3]*T+c.ns_dict[INP_type][4]*dt)
                # ??
             elif c.ns_dict[INP_type][0] == 2:
                 NS1[i] = 10e0**(c.ns_dict[INP_type][1]*T + c.ns_dict[INP_type][2]) 
            
             elif c.ns_dict[INP_type][0] == 3:
                 NS1[i] = 1e4*np.exp(c.ns_dict[INP_type][1]+c.ns_dict[INP_type][2]*T
                                 +c.ns_dict[INP_type][3]*T**2
                                 +c.ns_dict[INP_type][4]*T**3)
                 
             ACT_SITE[mode,:] += (np.pi *(6*M[mode,:]*Mass_frac[INP_type][mode]/
                                 (np.pi*c.aerosol_dict[INP_type][1]))**(2/3))*NS1[i]    
             #print(M[mode,:]*n.Mass_frac[INP_type][mode])
    NS_out = np.reshape(ACT_SITE,nbins*nmodes)      
    
    return NS_out
   
def icenucleation(MWAT, MBIN2, n_aer_bin, T, P, 
                  nbins, nmodes, rhobin, KAPPA, 
                  ncomps, dt, MBIN2_ICE, YICE_old, 
                  IND1, IND2, mass_for_act, RH, alpha_crit,
                  Heterogeneous_freezing_criteria, Mass_frac):
    """ calculate number of ice crystals nucleated this time step
        homogeneous aw -> Koop et al (2000) 
        heterogeneous ns -> Connolly et al (2012) """
        
    # calculate water activity
    RHOAT = MWAT/c.rhow+(MBIN2/rhobin)
    RHOAT = (MWAT+(MBIN2))/RHOAT
    Dw = ((MWAT + (MBIN2))*6/(np.pi*RHOAT))**(1/3)
    Dd = ((MBIN2*6)/(rhobin*np.pi))**(1/3)
    
    AW = (Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))
    
    # calculate homogeneous freezing rate following Koop et al (2000)
    # JW units m^-3 s^-1
    JW = KOOPNUCLEATIONRATE(AW, T, P, nbins, nmodes)
    
    # find the number of homogeneously frozen drops 
    #P=1-exp(-J*V*t), right hand column page 3 Koop et al
    number_frozen = np.absolute(-1*n_aer_bin * (np.exp(-JW * (MWAT/c.rhow) * dt)-1))
    #number_frozen = 0.0# switch off homogeneous freezing
    # update ice number
    YICE = YICE_old[IND1:IND2] + number_frozen +1e-50
    
    # update aerosol number
    n_aer_bin = n_aer_bin - number_frozen
    
    # aerosol mass going into bins
    DMAER = MBIN2_ICE*YICE_old[IND1:IND2]+YICE*MBIN2
    
    # mass of water going into bins
    M01 = YICE_old[IND1:IND2]*YICE_old[:IND1]+YICE*MWAT+1e-50
    
    # average MWAT of frozen particles
    M01 = np.where(M01 > 1e-50, M01/YICE, 1e-50) 
    
    # average aerosol mass in each bin
    MBIN2_ICE = DMAER/(YICE)
    
    # heteogeneous freezing
    NS = ACTIVESITES(T, MBIN2, rhobin, nbins, nmodes, ncomps,dt,Mass_frac)
    
    # heterogeneous freezing criteria
    MWAT_CRIT = alpha_crit/6*np.pi*Dd**3*c.rhow

    if Heterogeneous_freezing_criteria.lower() == 'activated drops':
        NS[:] = np.where(MWAT < mass_for_act , 0.0, NS) # this is not the problem for ice number agreeing with ACPIM
    elif Heterogeneous_freezing_criteria.lower() ==  'rh>1  [default]':
        if RH < 1.0:
            NS[:] = 0.0
    elif Heterogeneous_freezing_criteria.lower() == 'threshold water mass':
        NS[:] = np.where(MWAT < MWAT_CRIT, 0.0, NS) 
    else:
        print('heterogeneous freezing criteria '+ Heterogeneous_freezing_criteria.lower() + 'unknown')
        return
      
    #DN01 = (n_aer_bin)*(1e0-np.exp(-NS))
    
    # find change in ice concentration
    DN01 = [max(aer*(1-np.exp(-ns))-old, 0.0) for aer, ns, old in zip(n_aer_bin, NS, YICE_old[IND1:IND2])]
    
    # move aerosol from liquid distribution to ice distribution
    n_aer_bin = n_aer_bin - DN01
    
    # add number of heterogeneous ice onto homogeneous ice
    YICE = YICE + DN01
    # aerosol mass going into bins
    DMAER = MBIN2_ICE*YICE_old[IND1:IND2]+YICE*MBIN2
    
    # mass of water going into bins
    M01_2 = YICE_old[IND1:IND2]*YICE_old[:IND1]+YICE*MWAT
    
    
    M01_2 = np.where(M01_2 > 1e-50, M01_2/YICE, 1e-50) 
    # average MWAT of frozen particles
    #M01_2 = M01_2/(YICE) # is this adding extra water mass?
    
    # average aerosol mass in each bin
    MBIN2_ICE = DMAER/(YICE)
         
    return M01_2, MBIN2_ICE, YICE, n_aer_bin, number_frozen, MWAT_CRIT

def CAPACITANCE01(MWAT, PHI):
    PHI = 1.0 # this makes pyACPIM agree with ACPIM
    RHOICE = c.rhoi
    VOL = np.where(MWAT > 0, MWAT/RHOICE, 1e-50)
    
    A = (3e0*VOL/(4e0*np.pi*PHI))**(1e0/3e0)
    C = A*PHI
       
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
    FV = 1.0 # changing these has a significant impact on total ice mass 
    FH = 1.0 # and RH, need to add ventilation02 function to calculate FV and FH
  
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
    
    GAMMAICE = INHERENTGROWTH(T) #this is working
    return  AR * np.exp((GAMMAICE - 1)/(GAMMAICE + 2)*np.log(((MICE - MICEOLD)/
                                       RHO_DEP)))
   # AR = AR * np.exp((GAMMAICE - 1)/(GAMMAICE + 2)*np.log((QV + (MICE - MICEOLD)/
    #                                   RHO_DEP)/QV))
###############################################################################
    
def find_act_mass(aer_mass,T, nbins, nmodes,rhobin, Kappa):
    """ function the find the mass of water required for each particle
        activate
    """
    def kk02(NW):
        """ Kappa Koehler theory, Petters and Kriedenwies (2007)
    """
        MWAT = NW*c.mw
        RH_ACT = 0
        mass_bin_centre = brac
   # T = output[idx,ITEMP]
        mult = -1.0
    
        RHOAT = MWAT/c.rhow+(brac/rhobin3)
        RHOAT = (MWAT+(brac))/RHOAT

        Dw = ((MWAT + (mass_bin_centre))*6/(np.pi*RHOAT))**(1/3)
        Dd = ((mass_bin_centre*6)/(rhobin3*np.pi))**(1/3)
        KAPPA = kappabin3
    
        sigma = surface_tension(T)
        RH_EQ = mult*((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                 np.exp((4*sigma*c.mw)/(c.R*T*c.rhow*Dw)))-RH_ACT
  
        return RH_EQ

    act_mass = np.zeros([nbins*nmodes])
    
    for i in range(nmodes*nbins):
        rhobin3 = rhobin[i]
        kappabin3 = Kappa[i]
        brac = aer_mass[i]

        act_mass[i] = minimize_scalar(kk02,bracket=(brac*0.1, brac*50), method='brent', tol=1e-30).x*c.mw

    return act_mass

def calc_CDP_concentration(liq, aer_mass, ice, aer_mass_ice, T, IND1, IND2, rhobin):
    """ function to find the concentration of particles in CDP size bins
    """
    CDP_bins = [3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0,
       16.0, 18.0, 20.0, 22.0, 24.0, 26.0, 28.0, 30.0, 32.0, 34.0, 36.0,
       38.0, 40.0, 42.0, 44.0, 46.0, 48.0, 50.0]
    
    # FSSP bins
    #CDP_bins = [5,7,9,11,13,15,17,19,21,23,25,27,29,31,22,35,38,41,44,47]


    CDP_bins = [x*1e-6 for x in CDP_bins]

    CDP_conc_liq = np.zeros([len(CDP_bins)])
    CDP_conc_ice = np.zeros([len(CDP_bins)])
    CDP_conc_total = np.zeros_like(CDP_conc_liq)
    
    mass_w = liq[:IND1]
    mass_ice = ice[:IND1]
    
    num_bin_liq = liq[IND1:IND2]
    num_bin_ice = ice[IND1:IND2]
    
    RHOAT = (mass_w/c.rhow) + (aer_mass/rhobin)
    RHOAT = (mass_w+(aer_mass))/RHOAT                            
    DIAM = ((mass_w + (aer_mass))*6/(np.pi*RHOAT))**(1/3)
    
    for i in DIAM:
        x = bisect.bisect(CDP_bins, i)
        if x == 0: continue # stop everything going in first bin, which wont be seen by CDP
        if x > 29: x = 29
        CDP_conc_liq[x] = num_bin_liq[x]
    
    if T < 273.15:    
        RHOAT = (mass_ice/c.rhoi) + (aer_mass_ice/rhobin)
        RHOAT = (mass_ice+(aer_mass_ice))/RHOAT                            
        DIAM_ICE = ((mass_ice + (aer_mass_ice))*6/(np.pi*RHOAT))**(1/3)
        
        for i in DIAM_ICE:
            x = bisect.bisect(CDP_bins, i)
            if x == 0: continue # stop everything going in first bin, which wont be seen by CDP
            if x > 29: x = 29
            CDP_conc_ice[x] = num_bin_ice[x]
            
    CDP_conc_total = CDP_conc_liq + CDP_conc_ice 
    
    
    return CDP_conc_liq, CDP_conc_ice, CDP_conc_total

    
   
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    