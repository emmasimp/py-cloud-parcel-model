#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 29 16:28:42 2017

@author: mbexkes3
"""
# Add simelements package to path so can import modules
import sys
sys.path.insert(0, "/home/mbexkes3/Desktop/simelements") 

from numpy import ones, zeros, pi, exp, log, reshape, sqrt, tile
from scipy.optimize import brentq
from scipy.special import erfc

import constants as c
from functions import surface_tension


def setup_grid(rhobin, kappabin, rkm, Dlow, nbins, nmodes, RH, sig, NAER, D_AER,T):
    rhobin2 = reshape(rhobin, [nmodes, nbins])
    kappabin2 = reshape(kappabin, [nmodes, nbins])
    
    ########### Set-up water mass grid ############################################
    MWATGRID = ones((nbins+1)) # mass of water bin edges
    MWATCENTRE = ones((nbins)) # mass of water bin centre
    
    mass_bin_edges = zeros([nmodes,nbins]) # mass of aerosol bin edges
    mass_bin_centre = zeros([nmodes,nbins]) #mass of aerosol bin centre
    
    D = zeros([nmodes,nbins+1]) # dry diameter of aerosol bin edges
    NBIN = zeros([nmodes,nbins]) # number of aerosol in each bin
    MWAT = zeros([nmodes,nbins]) # mass of water each bin
    
    # calculate water mass at bin edges
    for i, dummy in enumerate(MWATGRID):
        MWATGRID[i] = pi/6*Dlow**3*c.rhow*2**(i/rkm)
    
    # calculate water mass at bin centers
    MWATCENTRE = 0.5*(MWATGRID[0:nbins]+MWATGRID[1:nbins+1])
    
    #------------------------------------------------------------------------------    
    def kk03(MBIN):    
        
         """function when root findind finds mass of aerosol(MBIN) corresponding
         to water mass (MWATGRID or MWATCENTRE)"""
         
         rhoat = dummy_vars/c.rhow + MBIN/rhobin
         rhoat = (dummy_vars+MBIN)/rhoat
         
         Dd = ((MBIN*6.)/(pi*rhobin))**(1./3.)
         Dw = (((dummy_vars+MBIN)*6.)/(pi*rhoat))**(1./3.)
     
         return mult*((Dw**3-Dd**3)/(Dw**3-Dd**3*(1-kappabin))*
                      exp((4*surface_tension(T)*c.mw)/(c.R*T*c.rhow*Dw)))-RH    
                      
      #------------------------------------------------------------------------------
    
    #values for root finding
    mult = 1
    
    # find aerosol mass that corresponds to water mass at bin edges and bin centres
    for k in range(nmodes):
        for i in range(nbins):
            dummy_vars = MWATGRID[i]
            rhobin = rhobin2[k,i]
            kappabin = kappabin2[k,i]
            #mass_bin_edges[k,i]=solveq.zbrent(kk03,1e-30,1e4, maxniter=1000, bracket=True)
            mass_bin_edges[k,i] = brentq(kk03,1e-30,1e4,xtol=1e-100, maxiter = 1000)#'dummy' in IO.f90 line 645
            dummy_vars = MWATCENTRE[i]
            MWAT[k,i] = dummy_vars
            mass_bin_centre[k,i]=brentq(kk03,1e-30,1e4,xtol=1e-100, maxiter = 1000)#'MBIN' in IO.f90 line 649
            D[k,i]=(6*mass_bin_edges[k,i]/(pi*rhobin))**(1/3)
    
        D[k,-1] = D[k,nbins-1]*100 # think this is to make sure there is a big bin to catch big particles
    
    #now put aerosol in bins
        NBIN[k,0:] = NAER[k]*0.5*(
             erfc(-(log(D[k,1:]/D_AER[k])/(sqrt(2)*sig[k])))-
             erfc(-(log(D[k,0:nbins]/D_AER[k])/(sqrt(2)*sig[k]))))
        #there is something wrong here - this is putting too much in first bin...
        NBIN[k,0]=NAER[k]*0.5*(
            (erfc(-(log(D[k,1]/D_AER[k]/(sqrt(2)*sig[k]))))))
        
    NBIN2 = reshape(NBIN,nmodes*nbins)
    MWATGRID = tile(MWATGRID, nmodes)
    MWAT = reshape(MWAT, nmodes*nbins)
    D = reshape(D, nbins*nmodes+nmodes)
    mass_bin_edges = reshape(mass_bin_edges, nmodes*nbins)
    mass_bin_centre = reshape(mass_bin_centre, nmodes*nbins)
    
    return NBIN2, MWATGRID, MWAT, D, mass_bin_edges, mass_bin_centre
