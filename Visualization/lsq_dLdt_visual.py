#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 16:13:03 2020
Least-squares regression to observed and simulated retreat rates
Based on JNB script lst_square_dLdt.py

@author: lizz
"""

import numpy as np
from scipy.stats import linregress
import pylab as plt


# Input file for observations
obs_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/observed_terminus_centroids.csv',delimiter=',',skiprows=1)
sim_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/dense_simulated_termini.csv',delimiter=',',skiprows=1)
cold_sim_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/dense_simulated_termini-min30Cice.csv',delimiter=',',skiprows=1)

# # For looping through the glacier ids
min_glacier_id = 1
max_glacier_id = 194


# Initialize arrays for dLdt and error in dLdt
dLdt_obs = []
dLdt_obs_err = []
dLdt_sim = []
dLdt_sim_err = []
dLdt_cold = []
dLdt_cold_err = []
bad_p = []

glacier_id_advancing = []
glacier_id_all = []
count = 0
for glacier_id in range(min_glacier_id,max_glacier_id+1):
    # Pull out data for each glacier id
    f1=obs_data[:,0]==glacier_id
    f2=sim_data[:,0]==glacier_id
    f3=cold_sim_data[:,0]==glacier_id
    # Only compute if we have more than 2 measurements
    if sum(f1)>2:
        # Least squares fit to observations
        m1,b,rvalue,pvalue1,err1=linregress(obs_data[f1,1],obs_data[f1,2])
        m2,b,rvalue,pvalue2,err2=linregress(sim_data[f2,1],sim_data[f2,2])
        m3,b,rvalue,pvalue3,err3=linregress(cold_sim_data[f3,1],cold_sim_data[f3,2])
        # if pvalue2>0.05:
        #     count +=1 # tally bad simulations
        #     continue
        if ~np.isnan(m1):
            dLdt_obs.append(m1);dLdt_obs_err.append(err1)
            dLdt_sim.append(m2);dLdt_sim_err.append(err2)
            dLdt_cold.append(m3);dLdt_cold_err.append(err3)
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



dLdt_range=np.linspace(-2,2,101)

dLdt_obs = np.array(dLdt_obs)
dLdt_sim = np.array(dLdt_sim)
dLdt_cold = np.array(dLdt_cold)
dLdt_obs_err = np.array(dLdt_obs_err)
dLdt_sim_err = np.array(dLdt_sim_err)
dLdt_cold_err = np.array(dLdt_cold_err)
bad_p = np.array(bad_p)

# Make scatter plot of sims vs obs with error bars
fig, ax = plt.subplots(1)
fig.clf()
# plt.plot(dLdt_range, dLdt_range, ls='-', color='Gainsboro', lw=2, alpha=0.4)
plt.fill_between(dLdt_range,dLdt_range,-2.1,color='Gainsboro',alpha=0.4)
# plt.errorbar(dLdt_obs,dLdt_sim,xerr=dLdt_obs_err,yerr=dLdt_sim_err,color='k',fmt='o',ecolor='lightgray',capsize=5)
# plt.errorbar(dLdt_obs[bad_p],dLdt_sim[bad_p],xerr=dLdt_obs_err[bad_p],yerr=dLdt_sim_err[bad_p],color='k',fmt='o',fillstyle='none', ecolor='lightgray',capsize=5)
# plt.errorbar(dLdt_obs[~bad_p],dLdt_sim[~bad_p],xerr=dLdt_obs_err[~bad_p],yerr=dLdt_sim_err[~bad_p],color='k',fmt='o', fillstyle='full', ecolor='lightgray',capsize=5)
plt.errorbar(dLdt_obs[bad_p],dLdt_cold[bad_p],xerr=dLdt_obs_err[bad_p],yerr=dLdt_cold_err[bad_p],color='dimgrey',fmt='o',fillstyle='none', ecolor='lightgrey',capsize=5)
plt.errorbar(dLdt_obs[~bad_p],dLdt_cold[~bad_p],xerr=dLdt_obs_err[~bad_p],yerr=dLdt_cold_err[~bad_p],color='dimgrey',fmt='o', fillstyle='full', ecolor='lightgrey',capsize=5)
plt.axhline(0, color='black',linestyle='--')
plt.axvline(0, color='black',linestyle='--')
plt.xlabel('Observed rate of length change (km/a)', fontsize=14)
plt.ylabel('Simulated rate of length change (km/a)', fontsize=14)
plt.text(-0.64,-1.75,'Overestimate retreat',weight='bold')
plt.text(-0.64,-0.2,'Underestimate retreat',weight='bold')
plt.axis([-0.65,0.05,-2,0.1])
plt.xticks([-0.6, -0.4, -0.2, 0])
plt.yticks([-1.8, -1.2, -0.6, 0])
plt.tick_params(axis='both', length=5, width=2, labelsize=12)
plt.show()
#plt.savefig('/Users/lizz/dLdt_compare.pdf',bbox_inches='tight')
