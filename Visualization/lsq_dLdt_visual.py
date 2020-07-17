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

# # Flotation condition comparison
# floatcondition_data = np.loadtxt('/Users/lizz/Desktop/flotation_simulated_termini.csv',delimiter=',',skiprows=1)
# outlier_gids,gls = floatcondition_data[:,0],floatcondition_data[:,1]

# # For looping through the glacier ids
min_glacier_id = 1
max_glacier_id = 194


# Initialize arrays for dLdt and error in dLdt
dLdt_obs = []
dLdt_obs_err = []
dLdt_sim = []
dLdt_sim_err = []
outlier = []
dLdt_float = {}
dLdt_float_err = {}

glacier_id_advancing = []
glacier_id_all = []
count = 0
for glacier_id in range(min_glacier_id,max_glacier_id+1):
# for glacier_id in glaciers_to_plot:
    # Pull out data for each glacier id
    f1=obs_data[:,0]==glacier_id
    f2=sim_data[:,0]==glacier_id
    # Only compute if we have more than 2 measurements??
    if sum(f1)>2:
        # Least squares fit to observations
        m1,b,rvalue,pvalue1,err1=linregress(obs_data[f1,1],obs_data[f1,2])
        m2,b,rvalue,pvalue2,err2=linregress(sim_data[f2,1],sim_data[f2,2])
        # This will filter out weird simulation data
        if pvalue2>0.05:
            print('Glacier ID',glacier_id,'p-value exceeds 0.15')
            count = count + 1
        else:
            if ~np.isnan(m1):
                dLdt_obs.append(m1);dLdt_obs_err.append(err1)
                dLdt_sim.append(m2);dLdt_sim_err.append(err2)
                # Stack glacier ids  of glaciers that pass the test
                glacier_id_all.append(glacier_id)

                # # Now compare with flotation condition
                # if glacier_id in outlier_gids:
                #     outlier.append(True)
                #     f3=floatcondition_data[:,0]==glacier_id
                #     m3,b,rvalue,pvalue3,err3=linregress(floatcondition_data[f3,1],floatcondition_data[f3,2])
                #     dLdt_float[glacier_id]=m3
                #     dLdt_float_err[glacier_id]=err3
                # else:
                #     outlier.append(False)
                #     continue



dLdt_range=np.linspace(-2,2,101)

dLdt_obs = np.array(dLdt_obs)
dLdt_sim = np.array(dLdt_sim)
dLdt_obs_err = np.array(dLdt_obs_err)
dLdt_sim_err = np.array(dLdt_sim_err)
outlier = np.array(outlier)
# fiddled=True means glaciers that have had a tongue removed
# fiddled=np.array(fiddled)

# Make scatter plot of sims vs obs with error bars
fig, ax = plt.subplots(1)
fig.clf()
# I think it makes it too busy distinguishing between tongues and no tongues
plt.fill_between(dLdt_range,dLdt_range,-2.1,color='Gainsboro',alpha=0.4)
plt.errorbar(dLdt_obs,dLdt_sim,xerr=dLdt_obs_err,yerr=dLdt_sim_err,color='k',fmt='o',ecolor='lightgray',capsize=5)
# plt.errorbar(dLdt_obs[outlier],dLdt_sim[outlier],xerr=dLdt_obs_err[outlier],yerr=dLdt_sim_err[outlier],color='k',fmt='o',fillstyle='none', ecolor='lightgray',capsize=5)
# plt.errorbar(dLdt_obs[~outlier],dLdt_sim[~outlier],xerr=dLdt_obs_err[~outlier],yerr=dLdt_sim_err[~outlier],color='k',fmt='o', fillstyle='full', ecolor='lightgray',capsize=5)
# plt.errorbar(dLdt_obs, dLdt_float, xerr=dLdt_obs_err, yerr=dLdt_float_err, color='b', fmt='.', ecolor='lightgray', capsize=5)
# for gid, x, y in zip(sorted(list(set(obs_data[:,0]))), dLdt_obs, dLdt_sim):
#     if gid in outlier_gids:
#         y2 = dLdt_float[gid]
        # plt.text(x, y, str(int(gid)), color='r', fontsize=12)
        # plt.plot((x, x), (y, y2), color='b')
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
