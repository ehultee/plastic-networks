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

# # Floating bits
# float_data = np.loadtxt('/Users/lizz/Desktop/flotation_lengths.csv',delimiter=',',skiprows=1)
# id,float_length = float_data[:,0],float_data[:,1]
# tongue_removed = id[float_length>1000.0] # glacier ID where a floating tongue was removed

# # For looping through the glacier ids
# min_glacier_id = 1
# max_glacier_id = 194
glaciers_to_plot=[2, 9, 36, 42, 66, 68, 69, 70, 72, 75, 76, 78, 98, 103, 105, 112, 116, 127, 137, 154, 164, 175]


# Initialize arrays for dLdt and error in dLdt
dLdt_obs = []
dLdt_obs_err = []
dLdt_sim = []
dLdt_sim_err = []

glacier_id_advancing = []
glacier_id_all = []
fiddled = []
count = 0
# for glacier_id in range(min_glacier_id,max_glacier_id+1):
for glacier_id in glaciers_to_plot:
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

                # # Now see if a floating tongue was removed
                # if glacier_id in tongue_removed:
                #     fiddled.append(True)
                # else:
                #     fiddled.append(False)



dLdt_range=np.linspace(-2,2,101)

dLdt_obs = np.array(dLdt_obs)
dLdt_sim = np.array(dLdt_sim)
dLdt_obs_err = np.array(dLdt_obs_err)
dLdt_sim_err = np.array(dLdt_sim_err)
# fiddled=True means glaciers that have had a tongue removed
fiddled=np.array(fiddled)

# Make scatter plot of sims vs obs with error bars
fig, ax = plt.subplots(1)
fig.clf()

# I think it makes it too busy distinguishing between tongues and no tongues
plt.fill_between(dLdt_range,dLdt_range,-2.1,color='Gainsboro',alpha=0.4)
plt.errorbar(dLdt_obs,dLdt_sim,xerr=dLdt_obs_err,yerr=dLdt_sim_err,color='k',fmt='.',ecolor='lightgray',capsize=5)
# plt.errorbar(dLdt_obs[fiddled],dLdt_sim[fiddled],xerr=dLdt_obs_err[fiddled],yerr=dLdt_sim_err[fiddled],color='k',fmt='.',ecolor='lightgray',capsize=5)
# plt.errorbar(dLdt_obs[~fiddled],dLdt_sim[~fiddled],xerr=dLdt_obs_err[~fiddled],yerr=dLdt_sim_err[~fiddled],color='r',fmt='.',markersize=3,ecolor='lightgray',capsize=5)
for gid, x, y in zip(sorted(list(set(obs_data[:,0]))), dLdt_obs, dLdt_sim):
    plt.text(x, y, str(gid), color='r', fontsize=12)
plt.axhline(0, color='black',linestyle='--')
plt.axvline(0, color='black',linestyle='--')
plt.xlabel('Observed rate of advance (km/a)')
plt.ylabel('Simulated rate of advance (km/a)')
plt.text(-0.64,-1.75,'Overestimate retreat',weight='bold')
plt.text(-0.64,-0.2,'Underestimate retreat',weight='bold')
plt.axis([-0.65,0.05,-2,0.1])
plt.xticks([-0.6, -0.4, -0.2, 0])
plt.yticks([-1.8, -1.2, -0.6, 0])
plt.show()
#plt.savefig('/Users/lizz/dLdt_compare.pdf',bbox_inches='tight')
