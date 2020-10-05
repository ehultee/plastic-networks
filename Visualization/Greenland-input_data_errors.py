#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Plot retreat rate differences by glacier ID to identify regional features
Created on Tue Sep 29 16:08:30 2020

@author: lizz
"""

import numpy as np
from scipy.stats import linregress
import pylab as plt


# Input file for observations
obs_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/observed_terminus_centroids.csv',delimiter=',',skiprows=1)
sim_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/dense_simulated_termini.csv',delimiter=',',skiprows=1)
bed_error_data = np.loadtxt('/Users/lizz/Documents/GitHub/Data_unsynced/Autonetwork_data_QC/All-bedtopo_errors.csv', delimiter=',', skiprows=1)
vel_error_data = np.loadtxt('/Users/lizz/Documents/GitHub/Data_unsynced/Autonetwork_data_QC/All-velocity_errors.csv', delimiter=',', skiprows=1)
flotation_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/flotation_lengths.csv',delimiter=',',skiprows=1)

# # For looping through the glacier ids
min_glacier_id = 1
max_glacier_id = 194


# Initialize arrays for dLdt and error in dLdt
dLdt_obs = []
dLdt_obs_err = []
dLdt_sim = []
dLdt_sim_err = []
bad_p = []
bed_err = []
vel_err = []
flotation_len = []

glacier_id_advancing = []
glacier_id_all = []
count = 0
for glacier_id in range(min_glacier_id,max_glacier_id+1):
    # Pull out data for each glacier id
    f1=obs_data[:,0]==glacier_id
    f2=sim_data[:,0]==glacier_id
    f3=bed_error_data[:,0]==glacier_id
    f4=vel_error_data[:,0]==glacier_id
    f5=flotation_data[:,0]==glacier_id
    # Only compute if we have more than 2 measurements
    if sum(f1)>2:
        # Least squares fit to observations
        m1,b,rvalue,pvalue1,err1=linregress(obs_data[f1,1],obs_data[f1,2])
        m2,b,rvalue,pvalue2,err2=linregress(sim_data[f2,1],sim_data[f2,2])
        # if pvalue2>0.05:
        #     count +=1 # tally bad simulations
        #     continue
        if ~np.isnan(m1):
            dLdt_obs.append(m1);dLdt_obs_err.append(err1)
            dLdt_sim.append(m2);dLdt_sim_err.append(err2)
            bed_err.append(bed_error_data[f3,2])
            vel_err.append(vel_error_data[f4,2])
            flotation_len.append(flotation_data[f5,1])
            # Stack glacier ids  of glaciers that pass the test
            glacier_id_all.append(glacier_id)
            if pvalue2>0.15: # Filter weird/oscillating simulations
                bad_p.append(True)
                count +=1
            # elif pvalue1>0.15: # Filter badly constrained obs rates
            #     bad_p.append(True)
            else:
                bad_p.append(False)
                continue


glacier_id_all=np.array(glacier_id_all)
dLdt_obs = np.array(dLdt_obs)
dLdt_sim = np.array(dLdt_sim)
dLdt_obs_err = np.array(dLdt_obs_err)
dLdt_sim_err = np.array(dLdt_sim_err)
bad_p = np.array(bad_p)
bed_err = np.array(bed_err)
vel_err = np.array(vel_err)
flotation_len = np.array(flotation_len)
dLdt_diff = dLdt_obs-dLdt_sim


fig1, ax1 = plt.subplots(1, figsize=(12,3))
ax1.axhline(y=0, ls='--', color='k')
ax1.axvline(x=90, ls='-', color='DarkGrey', alpha=0.3)
ax1.axvline(x=120, ls='-', color='DarkGrey', alpha=0.3)
ax1.axvline(x=176, ls='-', color='DarkGrey', alpha=0.3)
ax1.scatter(glacier_id_all, dLdt_diff, s=0.5*flotation_len, marker='p', facecolors='Grey', alpha=0.3, edgecolors='k', label='Removed floating tongue len.')
ax1.scatter(glacier_id_all, dLdt_diff, marker='d', s=bed_err, facecolors='Grey', alpha=0.4, edgecolors='k', label='Bed topo error')
ax1.scatter(glacier_id_all, dLdt_diff, s=10*vel_err, marker='s', facecolors='Grey', alpha=0.5, edgecolors='k', label='Velocity error')
ax1.scatter(glacier_id_all[~bad_p], dLdt_diff[~bad_p], marker='o', c='k')
ax1.scatter(glacier_id_all[bad_p], dLdt_diff[bad_p], marker='o', facecolors='none', edgecolors='k')
lgnd = ax1.legend()
for i in range(len(lgnd.legendHandles)): 
    lgnd.legendHandles[i]._sizes = [40]
ax1.set(xlabel='MEaSUREs Glacier ID', ylabel='dL/dt_obs - dL/dt_sim [km/a]', yticks=[-1, 0, 1, 2])
plt.tight_layout()
plt.show()