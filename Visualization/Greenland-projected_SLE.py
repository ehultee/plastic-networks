## Generate plots of sea level equivalent ice flux from 2006-2100 Greenland projections
## 5 Jun 2020  EHU
## Based on Greenland-hindcast-SLE.py

import numpy as np
import matplotlib.pyplot as plt
import glob
#from matplotlib.colors import LogNorm
from matplotlib import cm
## Special import for SERMeQ modules
import sys
sys.path.insert(0, 'Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *


##-------------------
### READING IN SIMS
##-------------------

testyears = arange(0, 100, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence', 
#'RCP4pt5', 
#'RCP8pt5'
)

#datemarker = '2019-02-08' #markers on filename to indicate date run
tempmarker = 'min10Cice' #and temperature of ice
timestepmarker = '99a_dt025a' #and total time and timestep

## Define which glaciers are in the simulated set
glacier_ids = range(1,195) #tell the function which MEaSUREs glacier IDs you want to process.
not_present = (93, 94, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 169) #glacier IDs missing from set
added_jan19 = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)
errors = (5, 18, 19, 29, 71, 92, 95, 97, 101, 107, 108, 120, 134) #glacier IDs that crashed in hindcasting 12 Mar 2019
errors_2020 = (17, 51, 100, 102, 106, 109, 110, 113, 115, 117, 118, 121, 141, 168, 171) #glacier IDs that showed network problems 21 May 2019
rmv = np.concatenate((not_present, errors, errors_2020))
for n in rmv:
    try:
        glacier_ids.remove(n)
    except ValueError:
        pass
glaciers_simulated = glacier_ids #to plot all
#glaciers_simulated = (3, 153, 175)

full_output_dicts = {}
for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for gid in glaciers_simulated:
        fn = glob.glob('Documents/GitHub/Data_unsynced/SERMeQ_output/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] #using glob * to select files of multiple run dates
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

def scenario_sorted_SLE(scenario_dictionary):
    sd = scenario_dictionary
    pernetwork_cumul_fx = []
    pernetwork_cumul_sle = []
    for j, gid in enumerate(glaciers_simulated):
        branch_fx = [np.nan_to_num(sd['GID{}'.format(gid)][k]['Terminus_flux']) for k in range(len(sd['GID{}'.format(gid)]))]
        total_fx = sum(branch_fx, axis=0)
        total_sle = (1E-12)*np.array(total_fx)/(361.8) #Gt ice/mm SLE conversion
        cumul_fx = np.cumsum(total_fx)
        cumul_sle = np.cumsum(total_sle)
        pernetwork_cumul_fx.append(cumul_fx)
        pernetwork_cumul_sle.append(cumul_sle)
        if max(cumul_sle)!=cumul_sle[-1]:
            print 'GID {}: SLE {} mm'.format(gid, cumul_sle[-1])
    p = np.asarray(pernetwork_cumul_sle)[np.asarray(pernetwork_cumul_sle)[:,-1].argsort()[::]]
    sorted_sle = np.cumsum(p, axis=0)
    return sorted_sle

projected_sle = scenario_sorted_SLE(full_output_dicts['persistence'])

##---------------------------
### PLOT THE ASSOCIATED FLUX
##---------------------------

## Settings for plots   
labels = [str(gid) for gid in glaciers_simulated] #set what the glaciers will be called in plotting.  Default is simply their MEaSUREs ID
markers = ['o', '.', ',', '^', 'd', '*']
cmap = cm.get_cmap('Greys')
#colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
colors = cmap(linspace(0.1, 0.9, num=len(glaciers_simulated)))
alt_colors = cm.get_cmap('winter')([0, 0.1, 0.3, 0.5, 0.7, 0.8])


plt.figure(figsize=(12,8))
for j in range(len(projected_sle)):
    color_idx = (np.abs(155*np.array([1, 1, 0.3, 0.5, 0.7, 0.8]) - j)).argmin() #replace first two selections as we have manually set them
    ms_selection = mod(gid, len(markers))
    plt.plot(testyears[::], projected_sle[j], linewidth=1, color=alt_colors[color_idx])
    plt.plot(testyears[::4], projected_sle[j][::4], linewidth=0, marker=markers[ms_selection], ms=10, color=alt_colors[color_idx])
    if j==0:
        plt.fill_between(testyears[::], y1=projected_sle[j], y2=0, color=alt_colors[color_idx], alpha=0.8)  
    else:
        plt.fill_between(testyears[::], y1=projected_sle[j], y2=projected_sle[j-1], color=alt_colors[color_idx], alpha=0.8)     
#plt.plot([0, 20, 40, 60], [0, 14, 28, 42], color='k', linewidth=1, ls='-', alpha=0.8) #GRACE linear trend
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xticks([0, 2, 4, 6, 8])
#plt.axes().set_xlim(0, 8.75)
#plt.axes().set_xticklabels(['2006', '2008', '2010', '2012', '2014'])
plt.show()