## Bar chart of retreat to encircle Greenland map
## 12 Aug 2019  EHU


from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
#import shapefile
import glob
#from matplotlib.colors import LogNorm
from matplotlib import cm
#from shapely.geometry import *
from scipy import interpolate
#from scipy.ndimage import gaussian_filter
#from scipy.stats import gaussian_kde
#from osgeo import gdal
#from mpl_toolkits.basemap import Basemap
import mpl_toolkits.basemap.pyproj as pyproj
#from plastic_utilities_v2 import *
from GL_model_tools import *
#from flowline_class_hierarchy import *


##-------------------
### READING IN SIMS
##-------------------


testyears = np.arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence', 
#'RCP4pt5', 
#'RCP8pt5'
)

#datemarker = '2019-02-08' #markers on filename to indicate date run
tempmarker = 'min10Cice' #and temperature of ice
timestepmarker = '8a_dt025a' #and total time and timestep

## Define which glaciers are in the simulated set
glacier_ids = range(1,195) #MEaSUREs glacier IDs to process.
not_present = (93, 94, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 169) #glacier IDs missing from set
added_jan19 = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)
errors = (5, 18, 19, 29, 71, 92, 95, 97, 100, 101, 102, 106, 107, 108, 110, 113, 117, 120, 121, 134, 168, 171) #glacier IDs that crashed in hindcasting 12 Mar 2019 *or* showed network problems 21 May 2019
rmv = np.concatenate((not_present, errors))
for n in rmv:
    try:
        glacier_ids.remove(n)
    except ValueError:
        pass
glaciers_simulated = glacier_ids #adjust this to reflect whether you want to examine the whole set or a subset


full_output_dicts = {}
for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for gid in glaciers_simulated:
        fn = glob.glob('Documents/GitHub/Data_unsynced/Hindcasted_networks/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] #using glob * to select files of multiple run dates
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

avg_sim_rates = []
for gid in glaciers_simulated:
    sim_termini = full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::]
    cumulative_dL = -1*float(sim_termini[-1]) #terminus positions given in arclength upglacier from initial terminus, so negative final terminus position corresponds to positive dL and vice versa
    cumulative_dt = float(testyears[-1])
    avg_dLdt = cumulative_dL/cumulative_dt
    avg_sim_rates.append(avg_dLdt)

sim_total = 0.00875*np.array(avg_sim_rates)


##-------------------
### PLOT BAR CHART
##-------------------

bin_edges = (0, -5, -10, -15)
binned_retreats = []
for s in sim_total:
    if s <= bin_edges[0] and s > bin_edges[1]:
        binned_retreats.append(1)
    elif s <= bin_edges[1] and s > bin_edges[2]:
        binned_retreats.append(2)
    elif s <= bin_edges[2] and s > bin_edges[3]:
        binned_retreats.append(3)
    elif s <= bin_edges[3]:
        binned_retreats.append(4)

breaks = (64, 100, 164, 200)
topleft = abs(np.array(glaciers_simulated) - breaks[0]).argmin()
topright = abs(np.array(glaciers_simulated) - breaks[1]).argmin()
lowright = abs(np.array(glaciers_simulated) - breaks[2]).argmin()
lowleft = abs(np.array(glaciers_simulated) - breaks[3]).argmin()

border_color = 'DarkGrey'
gid_tickspacing = 10

## Left border
plt.figure(figsize=(0.5, 6))
plt.barh(glaciers_simulated[0:topleft], width=binned_retreats[0:topleft], color=border_color) # go up to first break point
plt.axes().set_xlim(4, 0) # use "outward normal" for retreat
plt.axes().set_xticks([4, 3, 2, 1, 0]) # x-ticks only on first border
plt.axes().set_xticklabels([]) #no tick labels (will have legend)
plt.axes().set_ylim(0, breaks[0])
plt.axes().yaxis.set_ticks_position('right')
plt.axes().set_yticks(np.arange(0, breaks[0], gid_tickspacing)) # apply consistent spacing
plt.axes().patch.set_visible(False)
plt.axes().spines['left'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().tick_params(direction='out', width=2, top='off', left='off')
plt.show()

## Top border
plt.figure(figsize=(3.3, 0.5))
plt.bar(glaciers_simulated[topleft+1:topright], height=binned_retreats[topleft+1:topright], color=border_color)
plt.axes().set_ylim(0, 4)
plt.axes().set_yticks([]) # no y-ticks
plt.axes().set_xlim(breaks[0], breaks[1])
plt.axes().set_xticks(np.arange(ceil(round(breaks[0])/gid_tickspacing)*gid_tickspacing, ceil(round(breaks[1])/gid_tickspacing)*gid_tickspacing, gid_tickspacing))
plt.axes().patch.set_visible(False)
plt.axes().spines['left'].set_visible(False)
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().tick_params(direction='out', width=2, top='off', left='off', right='off')
plt.show()

## Right border
plt.figure(figsize=(0.5, 6))
plt.barh(glaciers_simulated[topright+1:lowright], width=binned_retreats[topright+1:lowright], color=border_color)
plt.axes().set_xlim(0, 4) # use "outward normal" for retreat
plt.axes().set_xticks([]) # no x-ticks
plt.axes().set_ylim(breaks[2], breaks[1]) #run from top downward
plt.axes().set_yticks(np.arange(ceil(round(breaks[1])/gid_tickspacing)*gid_tickspacing, ceil(round(breaks[2])/gid_tickspacing)*gid_tickspacing, gid_tickspacing))
plt.axes().patch.set_visible(False)
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().spines['bottom'].set_visible(False)
plt.axes().tick_params(direction='out', width=2, top='off', bottom='off', right='off')
plt.show()

## Bottom border
plt.figure(figsize=(3.3, 0.5))
plt.bar(glaciers_simulated[lowright+1:lowleft], height=binned_retreats[lowright+1:lowleft], color=border_color)
plt.axes().set_ylim(0, 4)
plt.axes().set_yticks([]) #no y-ticks
plt.axes().set_xlim(breaks[3], breaks[2])
plt.axes().set_xticks(np.arange(ceil(round(breaks[2])/gid_tickspacing)*gid_tickspacing, ceil(round(breaks[3])/gid_tickspacing)*gid_tickspacing, gid_tickspacing))
plt.axes().patch.set_visible(False)
plt.axes().spines['right'].set_visible(False)
plt.axes().spines['left'].set_visible(False)
plt.axes().spines['top'].set_visible(False)
plt.axes().tick_params(direction='out', width=2, top='off', left='off', right='off')
plt.show()

## Legend
plt.figure(figsize=(0.5, 1), frameon=False)
plt.barh((1, 3, 5, 7), width=(4, 3, 2, 1), color=border_color)
plt.axes().set_xticks([])
plt.axes().set_yticks([])
#plt.axes().patch.set_visible(False)
plt.axes().axis('off')
plt.show()