#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 16 10:59:23 2018

@author: mbexkes3
"""
from numpy import zeros, pi, exp, log, reshape, sqrt
from scipy.optimize import brentq, brent
from scipy.special import erfc, erfinv

import constants as c
import namelist as n
from functions import surface_tension


def setup_grid(rhobin, kappabin, rkm, Dlow, nbins, nmodes, RH, sig, NAER, D_AER,T):
  #  rhobin = n.rhoa
  #  kappabin = n.k
  #  nbins = n.nbins
  #  nmodes = n.nmodes
    
    RH = n.RHint
    sig = n.sig
    NAER = n.NAER
    D_AER = n.D_AER
    
    #set - up some variables
    D = zeros(n.nbins+1)
    NBIN = zeros([n.nmodes, n.nbins])
    MBIN = zeros([n.nmodes, n.nbins])
    MWAT = zeros([n.nmodes, n.nbins])


    def FIND_ERFCINV(x):
        return DUMMY_VARS-erfc(x)

    def kk02(MWAT2):
        """ Kappa Koehler theory, Petters and Kriedenwies (2007)
            
        """
        mass_bin_centre = MBIN[N_SEL, N_SEL2]
       
        rhobin3 = n.rhoa[N_SEL]
        kappabin3 = n.k[N_SEL]
        
       # T = 280
        
        RHOAT = MWAT2/c.rhow+(mass_bin_centre/rhobin3)
        RHOAT = (MWAT2+(mass_bin_centre))/RHOAT
        #RHOAT = c.rhow
        Dw = ((MWAT2 + (mass_bin_centre))*6/(pi*RHOAT))**(1/3)
        Dd = ((mass_bin_centre*6)/(rhobin3*pi))**(1/3)
        KAPPA = (mass_bin_centre/rhobin3*kappabin3)/(mass_bin_centre/rhobin3)
        
        sigma = surface_tension(T)
        RH_EQ = (((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-KAPPA))*
                     exp((4*sigma*c.mw)/(c.R*T*c.rhow*Dw)))-RH_ACT)
    
        return RH_EQ

    
    for I1 in range(nmodes):
        DUMMY = NAER[I1]/nbins
        D[0] = 5e-9
        NBIN[I1, 0] = NAER[I1]*0.5*erfc(-(log(D[0]/D_AER[I1])/sqrt(2)/sig[I1]))
        
        for J1 in range(1, nbins):
            D[J1] = (DUMMY+NAER[I1]*0.5*
                     erfc(-(log(D[J1-1]/D_AER[I1])/sqrt(2)/sig[I1])))
            DUMMY_VARS = 2/NAER[I1]*D[J1]
            DUMMY_VARS = min(DUMMY_VARS, 2-1e-10)
            DUMMY2 = brentq(FIND_ERFCINV, -1e10, 1e10, xtol=1e-30)
            D[J1] = exp(-DUMMY2*sig[I1]*sqrt(2)+log(D_AER[I1]))
            
            NBIN[I1, J1] = (NAER[I1]*(0.5*erfc(-log(D[J1]/D_AER[I1])/sqrt(2.0)/sig[I1])-
                0.5*erfc(-log(D[J1-1]/D_AER[I1])/sqrt(2.0)/sig[I1])))
        
            MBIN[I1,:] = pi/6*D[:-1]**3*n.rhoa[I1]
   
    for I1 in range(nmodes):
        for J1 in range(nbins):
            N_SEL = I1
            N_SEL2 = J1
          
            MULT = -1
         
            MULT = 1
            RH_ACT = min(RH, 0.999)
            MWAT[I1, J1] = brentq(kk02, -1e-30, 1e10, xtol=1e-30, maxiter=500)
            
    NBIN = reshape(NBIN, nmodes*nbins)
    MWAT = reshape(MWAT, nmodes*nbins)

    MBIN = reshape(MBIN, nmodes*nbins)
    return NBIN, MWAT, MBIN        