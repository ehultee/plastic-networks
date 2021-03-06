## Load-in and plotting of forward projections for automatically-selected Greenland networks 
## Generalization of Greenland-summary_plotting.py, developed Apr-Jun 2018 by EHU

#from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
#import collections
import glob
#from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib.patches as mpatches
#from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
## Special import for SERMeQ modules
import sys
sys.path.insert(0, 'Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *


###--------------------------------------
#### GLACIERS & SCENARIOS TO BE COMPARED
###--------------------------------------

glacier_ids = range(1,195) #tell the function which MEaSUREs glacier IDs you want to process.
not_present = (93, 94, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 169) #glacier IDs missing from set
added_jan19 = (139, 140, 141, 142, 143, 159, 161, 172, 173, 177)
errors = (5, 18, 19, 29, 71, 92, 95, 97, 101, 107, 108, 120, 134) #glacier IDs that crashed in hindcasting 12 Mar 2019
rmv = np.concatenate((not_present, errors))
for n in rmv:
    try:
        glacier_ids.remove(n)
    except ValueError:
        pass
#glaciers_simulated = glacier_ids #to plot all
glaciers_simulated = (2, 50, 60, 70, 80, 90, 100, 154, 155, 174, 176, 179)

testyears = arange(0, 100, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence', 
#'RCP4pt5', 
#'RCP8pt5'
)

tempmarker = 'min10Cice' #and temperature of ice
timestepmarker = '99a_dt025a' #and total time and timestep

full_output_dicts = {}

for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for gid in glaciers_simulated:
        fn = glob.glob('Documents/GitHub/Data_unsynced/SERMeQ_output/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] #using glob * to select files of multiple run dates
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

perscenario_SLE = compare_scenario_sle(full_output_dicts)

###--------------------------------------
#### PLOTTING
###--------------------------------------   
  
## Settings for plots   
labels = [str(gid) for gid in glaciers_simulated] #set what the glaciers will be called in plotting.  Default is simply their MEaSUREs ID
markers = ['o', '.', ',', '^', 'd', '*']
styles = ['-', ':', '-.', '-', '-', '-']
cmap = cm.get_cmap('winter')
scenario_colors = cm.get_cmap('Blues')([0.1, 0.3, 0.5, 0.7, 0.9])
#colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
colors = cmap(linspace(0.1, 0.9, num=len(glaciers_simulated)))
alt_colors = cm.get_cmap('Greys')([0.2, 0.3, 0.5, 0.7, 0.9])


####-------------------
##### INDIVIDUAL SCENARIO BREAKDOWNS
####-------------------

####terminus
plt.figure()
for j, gid in enumerate(glaciers_simulated):
    print gid
    term_positions = full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::]
    ms_selection = mod(j, len(styles))
    if max(term_positions)>20000:
        plt.plot(testyears, -0.001*np.array(term_positions), linewidth=2, color='r', linestyle=styles[ms_selection], label='{}'.format(gid))
    else:
        plt.plot(testyears, -0.001*np.array(term_positions), linewidth=2, color='Gainsboro', linestyle=styles[ms_selection])
    #plt.plot(testyears[::4], -0.001*np.array(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::])[::4], linewidth=0, marker=markers[ms_selection], ms=10, color=colors[j])
plt.legend(loc='lower left')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Terminus change [km]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_ylim(-100, 1)
#plt.axes().set_yticks([-75, -50, -25, 0])
plt.title('Terminus retreat of {} Greenland outlet glaciers 2006-2100 ERA-I, Tice=-10 C'.format(len(glaciers_simulated)), fontsize=20)
plt.show()
#
##
#####SINGLE NETWORK - splitting termini
#single_network_output = full_output_dicts['persistence']['GID10']
#plt.figure('Single network terminus change')
#for k in range(len(single_network_output)): #for each branch j
#    #colork = alt_colors[mod(k, len(colors))]
#    branch_termini = single_network_output[k]['Termini'] #for some reason enumerate doesn't work with loaded-in output, so we're stuck with this
#    #markerk = (k+2, mod(k+1, 3), 0)
#    plt.plot(testyears, -0.001*np.array(branch_termini[0:-1:]), linewidth=4, label='Branch {}'.format(k))
#    #plt.plot(testyears[::10], -0.001*np.array(branch_termini[0:-1:])[::10], linewidth=0, color=colork, marker=markerk, ms=10)
#plt.legend(loc='lower left')
#plt.axes().set_xlabel('Year of simulation', size=30)
#plt.axes().set_ylabel('Terminus change [km]', size=30)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_ylim(-50, 1)
##plt.axes().set_yticks([-50, -25, 0])
#plt.show()
#
#####Flux
#for j, gid in enumerate(glaciers_simulated):
#    plt.figure('GID{}'.format(gid))
#    ms_selection = mod(gid, len(styles))
#    plt.plot(testyears, 1E-12*np.array(pernetwork_cumul_fx[j]), linewidth=4, label=str(gid), color=colors[ms_selection], linestyle = styles[ms_selection])
#    plt.plot(testyears[::20], 1E-12*np.array(pernetwork_cumul_fx[j][::20]), linewidth=0, ms=10, marker=markers[ms_selection])
#    plt.fill_between(testyears, y1=1E-12*np.array(pernetwork_cumul_fx[j]), y2=0, alpha=0.5)
#    plt.axes().set_xlabel('Year of simulation', size=20)
#    plt.axes().set_ylabel('Cumulative ice flux [Gt]', size=20)
#    plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#    plt.axes().set_xlim(0,50)
#    plt.axes().set_xticks([0,10,20, 30, 40, 50])
#    #plt.axes().set_ylim(0, 401)
#    #plt.axes().set_yticks([0, 50, 100, 150, 200, 250, 300, 350, 400])
#    plt.show()
    ##
#######Sea level equivalent
plt.figure(figsize=(12,8))
for j, gid in enumerate(glaciers_simulated):
    #if gid!=10:
    ms_selection = mod(gid, len(styles))
    plt.plot(testyears[::], scenario_sle[j], linewidth=4, color=colors[ms_selection], label=gid)
    plt.plot(testyears[::5], scenario_sle[j][::5], linewidth=0, marker=markers[ms_selection], ms=10, color=colors[ms_selection])
    if j==0:
        plt.fill_between(testyears[::], y1=scenario_sle[j], y2=0, color=colors[ms_selection], alpha=0.7)  
    else:
        plt.fill_between(testyears[::], y1=scenario_sle[j], y2=scenario_sle[j-1], color=colors[ms_selection], alpha=0.7)     
#plt.plot([0, 20, 40, 60], [0, 14, 28, 42], color='k', linewidth=1, ls='-', alpha=0.8) #GRACE linear trend
#rect = mpatches.Rectangle((98,8.5), width=2, height=4.6, color='k', alpha=0.7) # Nick et al total projection for 2100, A1B
#rect2 = mpatches.Rectangle((98,11.3), width=2, height=6.2, color='k', alpha=0.7) # Nick et al total projection for 2100, RCP8.5 
#plt.axes().add_patch(rect)
#plt.axes().add_patch(rect2)
plt.legend(loc='upper left')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0, 100)
#plt.axes().set_xticks([0, 25, 50, 75, 100])
#plt.axes().set_ylim(0, 12)
#plt.axes().set_yticks([0, 2, 4, 6, 8, 10, 12])
plt.show()