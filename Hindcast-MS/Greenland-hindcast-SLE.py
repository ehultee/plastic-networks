## Generate plots of sea level equivalent ice flux from 2006-2014 Greenland hindcasts
## 13 Nov 2019  EHU
## Separate code from MEaSUREs-validation.py because *this cannot be directly validated*

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

testyears = arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
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
errors = (5, 17, 18, 19, 29, 51, 71, 92, 95, 97, 100, 101, 102, 106, 107, 108, 109, 110, 113, 115, 117, 120, 121, 134, 168, 171) #glacier IDs that crashed in hindcasting 12 Mar 2019 *or* showed network problems 21 May 2019
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


hindcast_sle = scenario_cumulative_SLE(full_output_dicts['persistence'])

##---------------------------
### PLOT THE ASSOCIATED FLUX
##---------------------------

## Settings for plots   
labels = [str(gid) for gid in glaciers_simulated] #set what the glaciers will be called in plotting.  Default is simply their MEaSUREs ID
markers = ['o', '.', ',', '^', 'd', '*']
styles = ['-', ':', '-.', '-', '-', '-']
cmap = cm.get_cmap('winter')
scenario_colors = cm.get_cmap('Blues')([0.1, 0.3, 0.5, 0.7, 0.9])
#colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
colors = cmap(linspace(0.1, 0.9, num=len(glaciers_simulated)))
alt_colors = cm.get_cmap('Greys')([0.2, 0.3, 0.5, 0.7, 0.9])


plt.figure(figsize=(12,8))
for j, gid in enumerate(glaciers_simulated):
    #if gid!=10:
    ms_selection = mod(gid, len(styles))
    plt.plot(testyears[::], hindcast_sle[j], linewidth=4, color=colors[ms_selection])
    plt.plot(testyears[::4], hindcast_sle[j][::4], linewidth=0, marker=markers[ms_selection], ms=10, color=colors[ms_selection])
    if j==0:
        plt.fill_between(testyears[::], y1=hindcast_sle[j], y2=0, color=colors[ms_selection], alpha=0.7)  
    else:
        plt.fill_between(testyears[::], y1=hindcast_sle[j], y2=hindcast_sle[j-1], color=colors[ms_selection], alpha=0.7)     
#plt.plot([0, 20, 40, 60], [0, 14, 28, 42], color='k', linewidth=1, ls='-', alpha=0.8) #GRACE linear trend
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
plt.axes().set_xticks([0, 2, 4, 6, 8])
plt.axes().set_xlim(0, 8.75)
plt.axes().set_xticklabels(['2006', '2008', '2010', '2012', '2014'])
plt.show()