## Load-in and plotting of forward projections for automatically-selected Greenland networks 
## Generalization of Greenland-summary_plotting.py, developed Apr-Jun 2018 by EHU

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import csv
import collections
#from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib.patches as mpatches
#from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from plastic_utilities_v2 import *
from GL_model_tools import *
from flowline_class_hierarchy import *



###-------------------
#### DEFINING NECESSARY FUNCTIONS
###-------------------
glacier_names = ('Sermeq Kujalleq [main]', 'Sermeq Kujalleq [central]', 'Sermeq Kujalleq [north]', 'Koge Bugt', 'Helheim', 'Kangerlussuaq')
ty_50a = arange(50, step=0.25)
ty_100a = arange(100, step=0.25)

### Load-in functionality to use only terminus position and flux
def lightload(filename, glacier_name, output_dictionary):
    output_dictionary[glacier_name] = {}
    
    with open(filename, 'rb') as handle:
        loadin = pickle.load(handle)
    
    N_Flowlines = loadin['N_Flowlines']
    mainline_termini = loadin['mainline_model_output']['Termini']
    mainline_flux = loadin['mainline_model_output']['Terminus_flux']
    output_dictionary[glacier_name][0] ={'Termini': mainline_termini, 'Terminus_flux': mainline_flux}
    
    if N_Flowlines >1:
        for n in range(N_Flowlines)[1::]:
            key_n = 'model_output_'+str(n)
            termini_n = loadin[key_n]['Termini']
            termflux_n = loadin[key_n]['Terminus_flux']
            output_dictionary[glacier_name][n] = {'Termini': termini_n, 'Terminus_flux': termflux_n}
    else:
        pass
        
    return output_dictionary

###--------------------------------------
#### GLACIERS & SCENARIOS TO BE COMPARED
###--------------------------------------

glaciers_simulated = (4, 8, 10, 11) #list MEaSUREs glacier IDs
testyears = arange(20, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence', 
#'RCP4pt5', 
#'RCP8pt5'
)

datemarker = '30Jan19' #markers on filename to indicate date run
tempmarker = 'min30Cice' #and temperature of ice
timestepmarker = '20a_dt025a' #and total time and timestep

full_output_dicts = {}

for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for i, gid in enumerate(glaciers_simulated):
        fn = 'GID{}-{}-{}-{}-{}.pickle'.format(gid, datemarker, s, tempmarker, timestepmarker)
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

perscenario_fluxes = []
perscenario_SLE = []

for s in full_output_dicts.keys():
    print 'Scenario {}'.format(s)
    out = full_output_dicts[s]
    pernetwork_cumul_fx = []
    pernetwork_cumul_sle = []
    for j, gid in enumerate(glaciers_simulated):
        branch_fx = [np.nan_to_num(out['GID{}'.format(gid)][k]['Terminus_flux']) for k in range(len(out['GID{}'.format(gid)]))]
        total_fx = sum(branch_fx, axis=0)
        total_sle = (1E-12)*np.array(total_fx)/(361.8) #Gt ice/mm SLE conversion
        cumul_fx = np.cumsum(total_fx)
        cumul_sle = np.cumsum(total_sle)
        pernetwork_cumul_fx.append(cumul_fx)
        pernetwork_cumul_sle.append(cumul_sle)
    scenario_flux = np.cumsum(pernetwork_cumul_fx, axis=0)
    perscenario_fluxes.append(scenario_flux[-1])
    scenario_sle = np.cumsum(pernetwork_cumul_sle, axis=0)
    print max(scenario_sle[-1])
    print(scenario_sle[-1][-1])
    perscenario_SLE.append(scenario_sle[-1])

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

plt.figure()
for j in range(len(perscenario_SLE)):
    plt.plot(testyears, perscenario_SLE[j], color=alt_colors[j+2], label=scenarios[j])
#for k in range(len(coarser_SLE)):
#    plt.plot(arange(100, step=0.5), coarser_SLE[k], color=scenario_colors[k+3], label=scenarios[k+3])
plt.legend(loc='upper left')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0, 100)
#plt.axes().set_xticks([0, 25, 50, 75, 100])
##plt.axes().set_ylim(-16, 1)
##plt.axes().set_yticks([0, 1, 2, 3, 4])
plt.show()

##### Scenario build
##for k in range(len(perscenario_SLE)):
##    plt.figure(k)
##    for j in range(k+1):
##        plt.plot(ty_100a, perscenario_SLE[j], color=alt_colors[j], label=scenarios[j], linewidth=4)
##    #for m in range(len(coarser_SLE)):
##    #    plt.plot(arange(100, step=0.5), coarser_SLE[m], color=alt_colors[m+3], label=scenarios[m+3])
##    plt.legend(loc='upper left')
##    plt.axes().set_xlabel('Year of simulation', size=20)
##    plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
##    plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##    plt.axes().set_xlim(0, 100)
##    plt.axes().set_xticks([0, 25, 50, 75, 100])
##    plt.axes().set_ylim(0, 40)
##    plt.axes().set_yticks([0, 5, 10, 15, 20, 25, 30, 35, 40])
##    plt.axes().set_yticklabels(['0', '', '10', '', '20', '', '30', '', '40'])
##    plt.show()
#
### Plotting warm in red and cool in blue
#plt.figure('Scenariodiff', figsize=(12, 8))
#plt.plot(ty_100a, perscenario_SLE[0], color='MediumBlue', label=scenarios[0], linewidth=2, ls='--') #cold baseline
#plt.plot(ty_100a, perscenario_SLE[1], color='FireBrick', label=scenarios[1], linewidth=2, ls='--') #warm baseline
#plt.plot(ty_100a, perscenario_SLE[2], color='MediumBlue', linewidth=2) #cold persistence
#plt.plot(ty_100a[::20], perscenario_SLE[2][::20], color='MediumBlue', label=scenarios[2], linewidth=0, marker='o', ms=6) #cold persistence
#plt.plot(ty_100a, perscenario_SLE[3], color='FireBrick', linewidth=2) #warm persistence
#plt.plot(ty_100a[::20], perscenario_SLE[3][::20], color='FireBrick', label=scenarios[3], linewidth=0, marker='o', ms=6) #warm persistence
#plt.plot(ty_100a, perscenario_SLE[4], color='FireBrick', linewidth=2) #warm doubling
#plt.plot(ty_100a[::20], perscenario_SLE[4][::20], color='FireBrick', label=scenarios[4], linewidth=0, marker='^', ms=8) #warm doubling
#plt.plot(ty_100a, perscenario_SLE[5], color='MediumBlue', linewidth=2) #cold persistence
#plt.plot(ty_100a[::20], perscenario_SLE[5][::20], color='MediumBlue', label=scenarios[5], linewidth=0, marker='.', ms=6) #cold persistence - modified to 0.8*forcing for SK main
#plt.plot([0, 20, 40, 60], [0, 14, 28, 42], color='k', linewidth=1, ls='-', alpha=0.8) #GRACE linear trend
#rect = mpatches.Rectangle((98,8.5), width=2, height=4.6, color='k', alpha=0.7) # Nick et al total projection for 2100, A1B
#rect2 = mpatches.Rectangle((98,11.3), width=2, height=6.2, color='k', alpha=0.7) # Nick et al total projection for 2100, RCP8.5
#plt.axes().add_patch(rect)
#plt.axes().add_patch(rect2)
#plt.legend(loc='upper left')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0, 100)
#plt.axes().set_xticks([0, 25, 50, 75, 100])
#plt.axes().set_ylim(0,30)
#plt.axes().set_yticks([0,10,20,30])
#plt.tight_layout()
#plt.show()
#
### Branch separation
##for j, s in enumerate(scenarios):
##    figname = 'Scenario_{}-branchsep'.format(s)
##    plt.figure(figname)
##    for k in range(len(output_dicts[j]['Helheim'])): #for each branch j
##        #colork = matplotlib.cm.get_cmap('viridis')(k/len(Helheim.model_output))
##        #markerk = (k+2, mod(k+1, 3), 0)
##        lsk = (':', '-.', '--', '-')
##        plt.plot(ty_100a[::], -0.001*np.array(CB_100a_output['Helheim'][k]['Termini'])[:-1:], linewidth=4, color='k', ls=lsk[k], label='Helheim line {}'.format(k))
##        #plt.plot(testyears[::10], -0.001*np.array(mo['Termini'])[:-1:10], linewidth=0, color='k', marker=markerk, ms=10, label='Helheim line {}'.format(k))
##    #plt.legend(loc='upper right')
##    plt.axes().set_xlabel('Year of simulation', size=20)
##    plt.axes().set_ylabel('Terminus change [km]', size=20)
##    plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##    plt.axes().set_xlim(0, 100)
##    plt.axes().set_xticks([0, 25, 50, 75, 100])
##    plt.axes().set_xticklabels(['0', '', '50', '', '100'])
##    #plt.axes().set_ylim(-16, 1)
##    plt.axes().set_yticks([-40, -30, -20, -10, 0])
##    plt.show()
##
##for j, s in enumerate(scenarios):
##    figname = 'Scenario_{}-branchflux'.format(s)
##    plt.figure(figname)
##    for k in range(len(output_dicts[j]['Helheim'])): #for each scenario j
##        #colork = matplotlib.cm.get_cmap('viridis')(k/len(Helheim.model_output))
##        #markerk = (k+2, mod(k+1, 3), 0)
##        lsk = (':', '-.', '--', '-')
##        plt.plot(ty_100a[::], np.cumsum(1E-12*np.nan_to_num(output_dicts[j]['Helheim'][k]['Terminus_flux']))[::], linewidth=4, color='k', ls=lsk[k], label='Helheim line {}'.format(k))
##        #plt.plot(testyears[::10], -0.001*np.array(mo['Termini'])[:-1:10], linewidth=0, color='k', marker=markerk, ms=10, label='Helheim line {}'.format(k))
##    plt.legend(loc='upper left')
##    plt.axes().set_xlabel('Year of simulation', size=20)
##    plt.axes().set_ylabel('Terminus flux [Gt]', size=20)
##    plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##    plt.axes().set_xlim(0, 100)
##    plt.axes().set_xticks([0, 25, 50, 75, 100])
##    plt.axes().set_xticklabels(['0', '', '50', '', '100'])
##    #plt.axes().set_ylim(-16, 1)
##    #plt.axes().set_yticks([-40, -30, -20, -10, 0])
##    plt.show()
##
##plt.figure('Summary')
##for j in range(len(perscenario_SLE)):
##    plt.plot(ty_100a, perscenario_SLE[j], color=alt_colors[j], label=scenarios[j], linewidth=4)
##for m in range(len(coarser_SLE)):
##    plt.plot(arange(100, step=0.5), coarser_SLE[m], color=alt_colors[m+3], label=scenarios[m+3], linewidth=4)
##plt.legend(loc='upper left')
##plt.axes().set_xlabel('Year of simulation', size=20)
##plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
##plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_xlim(0, 100)
##plt.axes().set_xticks([0, 25, 50, 75, 100])
##plt.axes().set_ylim(0, 35)
###plt.axes().set_yticks([0, 1, 2, 3, 4])
##plt.show()
#
#
#
#
####-------------------
##### INDIVIDUAL SCENARIO BREAKDOWNS
####-------------------
#
##Warm baseline fluxes
##Kanger_multibranch_flux = [Kanger_warmbaseline.model_output[j]['Terminus_flux'] for j in range(len(Kanger.flowlines))]
##Kanger_total_flux = sum(Kanger_multibranch_flux, axis = 0) #note that Kanger_multibranch_flux is multidimensional, needs care in summing
##Helheim_multibranch_flux = [Helheim_warmbaseline.model_output[j]['Terminus_flux'] for j in range(len(Helheim.flowlines))]
##Helheim_total_flux = sum(Helheim_multibranch_flux, axis=0)
##Jakobshavn_multibranch_flux = [Jak_main_warmbaseline.model_output[0]['Terminus_flux'], Jak_sec_warmbaseline.model_output[0]['Terminus_flux'], Jak_tert_warmbaseline.model_output[0]['Terminus_flux']]
##Jakobshavn_total_flux = sum(Jakobshavn_multibranch_flux, axis=0)
##fluxes = [Jak_main_warmbaseline.model_output[0]['Terminus_flux'], Jak_sec_warmbaseline.model_output[0]['Terminus_flux'], Jak_tert_warmbaseline.model_output[0]['Terminus_flux'], KogeBugt_warmbaseline.model_output[0]['Terminus_flux'], Helheim_warmbaseline.model_output[0]['Terminus_flux'], Kanger_warmbaseline.model_output[0]['Terminus_flux']]
##total_fluxes = [Jakobshavn_total_flux, KogeBugt_warmbaseline.model_output[0]['Terminus_flux'], Helheim_total_flux, Kanger_total_flux]
##NOTE NEED TO INCLUDE FLUXES FROM ALL LINES OF HELHEIM, KANGER
#
#
####Cold baseline fluxes for plot
##Kanger_multibranch_flux_cold = [np.nan_to_num(CB_100a_output['Kangerlussuaq'][j]['Terminus_flux']) for j in range(len(CB_100a_output['Kangerlussuaq']))]
##Kanger_total_flux_cold = sum(Kanger_multibranch_flux_cold, axis = 0) #note that Kanger_multibranch_flux is multidimensional, needs care in summing
##Helheim_multibranch_flux_cold = [np.absolute(CB_100a_output['Helheim'][j]['Terminus_flux']) for j in range(len(CB_100a_output['Helheim']))]
##Helheim_total_flux_cold = sum(Helheim_multibranch_flux_cold, axis=0)
##total_fluxes_split_CB = [CB_100a_output['Sermeq Kujalleq [main]'][0]['Terminus_flux'], CB_100a_output['Sermeq Kujalleq [central]'][0]['Terminus_flux'], CB_100a_output['Sermeq Kujalleq [north]'][0]['Terminus_flux'], CB_100a_output['Koge Bugt'][0]['Terminus_flux'], Helheim_total_flux_cold, Kanger_total_flux_cold] #totals from all termini, with Jak split by branch (network)
##
##fluxes_cleaned = []
##sle = [] #will be array of annual sea level contributions
##for flux in total_fluxes_split_CB:
##    flux_c = np.nan_to_num(flux)
##    #flux_abs = np.absolute(flux_c) #quick fix 9 Jun for some lines that have negative width upstream.  Fix more substantially later.
##    fluxes_cleaned.append(flux_c)
##    sleq = (1E-12)*np.array(flux_c)/(361.8) #Gt ice / mm sea level equiv conversion
##    sle.append(sleq)
##cumul_sle_pernetwork = []
##total_sle = []
##for sl in sle:
##    c = np.cumsum(sl)
##    cumul_sle_pernetwork.append(c)
##total_sle = np.cumsum(cumul_sle_pernetwork, axis=0)
#
### Cold persistence fluxes for plot
#Kanger_multibranch_flux_CP = [np.nan_to_num(CP_100a_output['Kangerlussuaq'][j]['Terminus_flux']) for j in range(len(CP_100a_output['Kangerlussuaq']))]
#Kanger_total_flux_CP = sum(Kanger_multibranch_flux_CP, axis = 0) #note that Kanger_multibranch_flux is multidimensional, needs care in summing
#Helheim_multibranch_flux_CP = [np.absolute(CP_100a_output['Helheim'][j]['Terminus_flux']) for j in range(len(CP_100a_output['Helheim']))]
#Helheim_total_flux_CP = sum(Helheim_multibranch_flux_CP, axis=0)
#total_fluxes_split_CP = [CP_100a_output['Sermeq Kujalleq [main]'][0]['Terminus_flux'], CP_100a_output['Sermeq Kujalleq [central]'][0]['Terminus_flux'], CP_100a_output['Sermeq Kujalleq [north]'][0]['Terminus_flux'], CP_100a_output['Koge Bugt'][0]['Terminus_flux'], Helheim_total_flux_CP, Kanger_total_flux_CP] #totals from all termini, with Jak split by branch (network)
#
#fluxes_cleaned = []
#sle = [] #will be array of annual sea level contributions
#for flux in total_fluxes_split_CP:
#    flux_c = np.nan_to_num(flux)
#    #flux_abs = np.absolute(flux_c) #quick fix 9 Jun for some lines that have negative width upstream.  Fix more substantially later.
#    fluxes_cleaned.append(flux_c)
#    sleq = (1E-12)*np.array(flux_c)/(361.8) #Gt ice / mm sea level equiv conversion
#    sle.append(sleq)
#cumul_sle_pernetwork = []
#total_sle = []
#for sl in sle:
#    c = np.cumsum(sl)
#    cumul_sle_pernetwork.append(c)
#total_sle = np.cumsum(cumul_sle_pernetwork, axis=0)
#
#
##projections_WB = [Jak_main_warmbaseline.model_output, Jak_sec_warmbaseline.model_output, Jak_tert_warmbaseline.model_output, KogeBugt_warmbaseline.model_output, Helheim_warmbaseline.model_output, Kanger_warmbaseline.model_output]
##projections_CB = [Jak_main_coldbaseline.model_output, Jak_sec_coldbaseline.model_output, Jak_tert_coldbaseline.model_output, KogeBugt_coldbaseline.model_output, Helheim_coldbaseline.model_output, Kanger_coldbaseline.model_output]
#
#
###terminus
plt.figure()
for j, gid in enumerate(glaciers_simulated):
    print gid
    plt.plot(testyears, -0.001*np.array(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::]), linewidth=4, color=colors[j], linestyle=styles[j], label='{}'.format(gid))
    plt.plot(testyears[::4], -0.001*np.array(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::])[::4], linewidth=0, marker=markers[j], ms=10, color=colors[j])
plt.legend(loc='lower left')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Terminus change [km]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_ylim(-100, 1)
#plt.axes().set_yticks([-75, -50, -25, 0])
plt.title('20 a (dt=0.25 a) climate persistence - Gaussian sigma=2 - A=3.7E-26', fontsize=20)
plt.show()
#
##
##SINGLE NETWORK - splitting termini
single_network_output = full_output_dicts['persistence']['GID9']
plt.figure('Single network terminus change')
for k in range(len(single_network_output)): #for each branch j
    #colork = alt_colors[mod(k, len(colors))]
    branch_termini = single_network_output[k]['Termini'] #for some reason enumerate doesn't work with loaded-in output, so we're stuck with this
    #markerk = (k+2, mod(k+1, 3), 0)
    plt.plot(testyears, -0.001*np.array(branch_termini[0:-1:]), linewidth=4, label='Branch {}'.format(k))
    #plt.plot(testyears[::10], -0.001*np.array(branch_termini[0:-1:])[::10], linewidth=0, color=colork, marker=markerk, ms=10)
plt.legend(loc='lower left')
plt.axes().set_xlabel('Year of simulation', size=30)
plt.axes().set_ylabel('Terminus change [km]', size=30)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_ylim(-50, 1)
#plt.axes().set_yticks([-50, -25, 0])
plt.show()
##
#####Flux
##plt.figure()
##for j in range(len(names)):
##    plt.plot(testyears, 1E-12*np.array(fluxes[j]), linewidth=4, linestyle=styles[j], color=colors[j], label=names[j])
##    plt.plot(testyears[::50], 1E-12*np.array(fluxes[j][::50]), linewidth=0, marker=markers[j], ms=10, color=colors[j])
##    plt.fill_between(testyears, y1=1E-12*np.array(fluxes[j]), y2=0, color=colors[j], alpha=0.5)    
##plt.legend(loc='upper right')
##plt.axes().set_xlabel('Year of simulation', size=20)
##plt.axes().set_ylabel('Terminus ice flux [Gt/a]', size=20)
##plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_xlim(0,50)
##plt.axes().set_xticks([0,10,20, 30, 40, 50])
##plt.axes().set_ylim(0, 401)
##plt.axes().set_yticks([0, 50, 100, 150, 200, 250, 300, 350, 400])
##plt.show()
##
######Sea level equivalent
#plt.figure(figsize=(12,8))
#for j in range(len(names)):
#    plt.plot(1+ty_100a[::], total_sle[j], linewidth=4, color=colors[::][j], label=names[j])
#    plt.plot(1+ty_100a[::5], total_sle[j][::5], linewidth=0, marker=markers[j], ms=10, color=colors[::][j])
#    if j==0:
#        plt.fill_between(1+ty_100a[::], y1=total_sle[j], y2=0, color=colors[::][j], alpha=0.7)  
#    else:
#        plt.fill_between(1+ty_100a[::], y1=total_sle[j], y2=total_sle[j-1], color=colors[::][j], alpha=0.7)     
#plt.plot([0, 20, 40, 60], [0, 14, 28, 42], color='k', linewidth=1, ls='-', alpha=0.8) #GRACE linear trend
#rect = mpatches.Rectangle((98,8.5), width=2, height=4.6, color='k', alpha=0.7) # Nick et al total projection for 2100, A1B
#rect2 = mpatches.Rectangle((98,11.3), width=2, height=6.2, color='k', alpha=0.7) # Nick et al total projection for 2100, RCP8.5 
#plt.axes().add_patch(rect)
#plt.axes().add_patch(rect2)
#plt.legend(loc='upper left')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_xlim(0, 100)
#plt.axes().set_xticks([0, 25, 50, 75, 100])
#plt.axes().set_ylim(0, 12)
#plt.axes().set_yticks([0, 2, 4, 6, 8, 10, 12])
#plt.show()