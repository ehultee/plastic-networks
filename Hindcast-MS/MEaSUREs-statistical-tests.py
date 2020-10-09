#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Statistical tests for SERMeQ
Created on Wed May 27 21:39:28 2020
Namespace of Greenland-annual_dL_comparison for now

@author: EHU
"""
import numpy as np
from scipy.stats import spearmanr, kendalltau

# Input file for obs/sim comparison
obs_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/observed_terminus_centroids.csv',delimiter=',',skiprows=1)
sim_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/simulated_termini.csv',delimiter=',',skiprows=1)
dense_sim_data = np.loadtxt('/Users/lizz/Desktop/Hindcasted_networks/dense_simulated_termini.csv',delimiter=',',skiprows=1)

obs_by_gid = {}
sim_by_gid = {}
for o,s in zip(obs_data, sim_data):
    if o[0] not in obs_by_gid.keys():
        k= o[0]
        obs_by_gid[k] = [o[2]]
        sim_by_gid[k] = [s[2]]
    else:
        obs_by_gid[k].append(o[2])
        sim_by_gid[k].append(s[2])
        
## Regular set--"un-futzed" glaciers
corrs = []
spearmans = []
kendalls = []
zero_compare = []
for gid in glaciers_to_plot:
    tc = -1000*termpos_corrections[gid]
    sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
    st = -0.001*(tc+np.array(sim_termini))
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    otc = -1*obs_term_centr
    corrs.append(np.correlate(st, otc))
    spearmans.append(spearmanr(st, otc)[0])
    k, p = kendalltau(st, otc)
    kendalls.append(k)
    zero_compare.append(spearmanr(np.zeros(len(otc)), otc)[0])
c = np.squeeze(corrs)
s = np.squeeze(spearmans)
ks = np.squeeze(kendalls)
zc = np.squeeze(zero_compare)


## Futzed glaciers - compare statistics
corrs_f = []
spearmans_f = []
for gid in (set(added_jan19) | set(seaward_projected)): #set union of fussy glaciers
    if gid not in rmv:
        tc = -1000*termpos_corrections[gid]
        sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
        st = -0.001*(tc+np.array(sim_termini))
        obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
        obs_term_centr = obs_termini[:,1]
        otc = -1*obs_term_centr
        corrs_f.append(np.correlate(st, otc))
        spearmans_f.append(spearmanr(st, otc)[0])
cf = np.squeeze(corrs_f)
sf = np.squeeze(spearmans_f)

## All together
spearmans_t = []
spearmans_p = []
kendalls_t = []
for gid in (set(glaciers_to_plot) | set(added_jan19) | set(seaward_projected)): #set union of fussy glaciers
    if gid not in rmv:
        tc = -1000*termpos_corrections[gid]
        sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
        st = -0.001*(tc+np.array(sim_termini))
        obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
        obs_term_centr = obs_termini[:,1]
        otc = -1*obs_term_centr
        s, p = spearmanr(st, otc)
        # spearmans_t.append(spearmanr(st, otc)[0])
        spearmans_t.append(s)
        spearmans_p.append(p)
        kendalls_t.append(kendalltau(st, otc)[0])
st = np.squeeze(spearmans_t)
sp = np.squeeze(spearmans_p)
kt = np.squeeze(kendalls_t)


fig, ax = plt.subplots(1)
_,_, bars = ax.hist(st[~np.isnan(st)], bins=[-1, -0.75, -0.5, -0.25, 0, 0.25, 0.5, 0.75, 1.0], color='k')
for bar in bars:
    if bar.get_x()<0:
        bar.set_facecolor('LightGrey')
        bar.set_hatch('/')
ax.set_xlabel(r'Rank correlation $\rho$, observed vs. simulated termini', fontsize=18)
ax.set_ylabel('Count', fontsize=18)
ax.set_xticks([-1, -0.5, 0, 0.5, 1])
ax.tick_params(axis='both', length=5, width=2, labelsize=16)
plt.show()