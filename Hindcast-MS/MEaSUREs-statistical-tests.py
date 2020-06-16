#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Statistical tests for SERMeQ
Created on Wed May 27 21:39:28 2020
Namespace of Greenland-annual_dL_comparison for now

@author: EHU
"""

# Input file for obs/sim comparison
obs_data = np.loadtxt('/Users/lizz/Desktop/observed_terminus_centroids.csv',delimiter=',',skiprows=1)
sim_data = np.loadtxt('/Users/lizz/Desktop/simulated_termini.csv',delimiter=',',skiprows=1)
dense_sim_data = np.loadtxt('/Users/lizz/Desktop/dense_simulated_termini.csv',delimiter=',',skiprows=1)

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
kendalls_t = []
for gid in (set(glaciers_to_plot) | set(added_jan19) | set(seaward_projected)): #set union of fussy glaciers
    if gid not in rmv:
        tc = -1000*termpos_corrections[gid]
        sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
        st = -0.001*(tc+np.array(sim_termini))
        obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
        obs_term_centr = obs_termini[:,1]
        otc = -1*obs_term_centr
        spearmans_t.append(spearmanr(st, otc)[0])
        kendalls_t.append(kendalltau(st, otc)[0])
st = np.squeeze(spearmans_t)
kt = np.squeeze(kendalls_t)


