#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 12 15:56:34 2018

@author: mbexkes3
"""
from numpy import zeros, zeros_like, reshape

import namelist as n
import constants as c
import setup_grids as s

import importlib # NOTE: need to be in the correct directory for this to work
importlib.reload(n)
importlib.reload(c)
nmodes = n.nmodes
nbins = n.nbins

# -------------------- create index values -----------------------------------
IND1 = nbins*nmodes # index for mass
IND2 = IND1*2 # index for number
IPRESS = -3
ITEMP = -2
IRH = -1

IPRESS_ICE = -3
ITEMP_ICE = -2
IRH_ICE = -1
# ----------------------------------------------------------------------------

# --------------------- initailizes arrays ------------------------------------ 
t_final = n.runtime
Y = zeros(nbins*nmodes*2+3)
Y_AER = zeros([int(t_final/n.dt),IND1]) # mass of aerosol bin
YICE_AER = zeros_like(Y_AER) # mass of aerosol in ice bins

Y[IPRESS] = n.P
Y[ITEMP] = n.T
Y[IRH] = n.RH

rhobin = zeros(nmodes*nbins)
kappabin = zeros(nmodes*nbins)   
molwbin = zeros(nmodes*nbins) 

k = zeros(nmodes)
rhoa = zeros(nmodes)
molw_aer = zeros(nmodes)

Kappa = zeros([nmodes,nbins])# use this one in all calculations except setup grids

output = zeros([int(t_final/n.dt),len(Y)])
output_ice = zeros([int(t_final/n.dt),IND1*3+3])

CDP_CONC_liq = zeros([int(t_final/n.dt),30])
CDP_CONC_ice = zeros([int(t_final/n.dt),30])
CDP_CONC_total = zeros([int(t_final/n.dt),30])

dummy  = zeros_like(Y)
dummy2 = zeros_like(Y)
dummy4 = zeros_like(Y)

ice_aer = zeros([int(t_final/n.dt),IND1])

ACT_DROPS = zeros([int(t_final/n.dt),IND1])
ncomps = len(n.Mass_frac)


# Add in new aerosol types from user ------------------------------------------
if ncomps > len(c.aerosol_dict):
    for key in n.Mass_frac:
        if key not in c.aerosol_dict:
            if sum(n.Mass_frac[key][:]) > 0.0:
                print('Aerosol type', key, 'not found')
                print('please enter values for new aerosol type ')
                new_aerosol_type = list(range(4))
                print('Molecular weight in kg per mole ')
                new_aerosol_type[0] = input()
                print('Density kg/m^3 ')
                new_aerosol_type[1] = input()
                print('Kappa')
                new_aerosol_type[3] = input()
                #create new temporary dictionary entry 
                c.aerosol_dict[key] = list(map(float, new_aerosol_type))
                # add IN properties to new aerosol type
                print('Is aerosol type',key,'an INP?')
               
                INP = input('yes / no   ')
                if INP[0] == 'y' : 
                    is_INP = True
                else:
                    is_INP = False
                if is_INP:
                    new_INP_type = zeros(3)
                    print('Please provide ns values for ',key)
                    new_INP_type[0] = 2 # just use feldspar type at the moment
                    print('ns value 1')
                    new_INP_type[1] = input()
                    print('ns value 2')
                    new_INP_type[2] = input()
                    c.ns_dict[key] = list(map(float, new_INP_type))
    
                print('do you want to save new aerosol type ',key,' for future use')
                to_save = input('yes / no   ')
                if to_save[0] == 'y': 
                    to_save = True
                else:
                    to_save = False
                if to_save:

                    with open("constants.py", "r") as in_file:
                        buf = in_file.readlines()
                 
                    with open("constants.py", "w") as out_file:
                       for line in buf:
                         if line == "aerosol_dict = {'ammonium sulphate': [132.14e-3, 1770, 3, 0.61],\n":
                             line = line + "\'"+key+"\': ["+str(new_aerosol_type[0])+','+str(new_aerosol_type[1])+','+str(new_aerosol_type[2])+','+str(new_aerosol_type[3])+"],"
                         if line == "ns_dict = {'Kaolinite':[0, -0.8881, 12.9445, -0.8802, 231.383],\n":
                             if is_INP:
                                 line = line + "\'"+key+"\': ["+str(new_INP_type[0])+','+str(new_INP_type[1])+','+str(new_INP_type[2])+"],"
                         
                         out_file.write(line)
    # -----------------------------------------------------------------------------
# --------- catch error if sum of mass fractions is greater than 1 ------------
total = zeros(n.nmodes)

for mode in range(n.nmodes):
    for key in n.Mass_frac.keys():
        total[mode] += n.Mass_frac[key][mode]
    if total[mode] != 1:
        print('ERROR - Simulation Stopped')
        print('sum of mass fractions in mode', mode+1, 'does not equal 1')
        ERROR_FLAG = True
        break
    else:
        ERROR_FLAG = False
        continue

        
# -----------------------------------------------------------------------------
 
 # calculate average values for each mode
for mode in range(nmodes):
    k[mode] = sum([n.Mass_frac[key][mode]*c.aerosol_dict[key][3] for key in c.aerosol_dict.keys()])
    rhoa[mode] = 1/sum([n.Mass_frac[key][mode]/c.aerosol_dict[key][1] for key in c.aerosol_dict.keys()])
    molw_aer[mode] = sum([n.Mass_frac[key][mode]*c.aerosol_dict[key][0] for key in c.aerosol_dict.keys()])
                                    
for i in range(nmodes):
    rhobin[i*nbins:nbins+nbins*i] = rhoa[i]
    kappabin[i*nbins:nbins+nbins*i] = k[i]
    molwbin[i*nbins:nbins+nbins*i] = molw_aer[i]
    
############################### set-up grid ###################################
GRID = s.setup_grid(rhobin, kappabin, n.rkm, n.Dlow, nbins, 
                    nmodes, n.RH, n.sig, n.NAER, n.D_AER, n.T, rhoa, k)

Y[nbins*nmodes:-3] = GRID[0] # number of aerosol in each bin
Y[:IND1] = GRID[1]               # water mass in each bin
Y_AER[0,:] = GRID[2]                  # mass of aerosol in each bin
YICE_AER[0,:] = GRID[2]                 # mass of aerosol in ice bins 

for mode in range(nmodes):
    ind1 = nbins*mode
    ind2 = ind1+nbins
    Kappa[mode] = sum([n.Mass_frac[key][mode]*Y_AER[0,ind1:ind2]/
                         c.aerosol_dict[key][1]*c.aerosol_dict[key][3] for key
                         in c.aerosol_dict.keys()])/sum([n.Mass_frac[key][mode]
                         *Y_AER[0,ind1:ind2]/c.aerosol_dict[key][1] for key
                         in c.aerosol_dict.keys()])

Kappa = reshape(Kappa,[nbins*nmodes])



        











