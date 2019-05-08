## Plotting individual glacier results for case studies
## 24 Apr 2019  EHU

import numpy as np
import matplotlib.pyplot as plt
import csv
import shapefile
#import collections
import glob
#from matplotlib.colors import LogNorm
from matplotlib import cm
import matplotlib.patches as mpatches
from shapely.geometry import *
from scipy import interpolate
from scipy.ndimage import gaussian_filter
from plastic_utilities_v2 import *
from GL_model_tools import *
from flowline_class_hierarchy import *



###-------------------
#### DEFINING NECESSARY FUNCTIONS
###-------------------

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
#### GLACIERS TO PLOT
###--------------------------------------

gids_by_name = ((3, 'Jakobshavn/SK'), 
(175, 'Helheim'), 
(153, 'Kangerlussuaq'),
(10, 'Store'),
) #array pairing glacier IDs with their recognized  names

glaciers_to_plot = [g[0] for g in gids_by_name if g[1] in ('Jakobshavn/SK', 'Helheim', 'Kangerlussuaq', 'Store')] # select gids of glaciers to plot as case studies

testyears = arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence', 
#'RCP4pt5', 
#'RCP8pt5'
)

#datemarker = '2019-02-08' #markers on filename to indicate date run
tempmarker = 'min10Cice' #and temperature of ice
timestepmarker = '8a_dt025a' #and total time and timestep

full_output_dicts = {}

for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for gid in glaciers_to_plot:
        fn = glob.glob('Documents/GitHub/Data_unsynced/Hindcasted_networks/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] #using glob * to select files of multiple run dates
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    #for i, gid in enumerate(glaciers_to_plot):
    #    fn = 'GID{}-{}-{}-{}-{}.pickle'.format(gid, datemarker, s, tempmarker, timestepmarker)
    #    lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

perscenario_fluxes = []
perscenario_SLE = []

for s in full_output_dicts.keys():
    print 'Scenario {}'.format(s)
    out = full_output_dicts[s]
    pernetwork_cumul_fx = []
    pernetwork_cumul_sle = []
    for j, gid in enumerate(glaciers_to_plot):
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

## Reading in observed termini and saved networks, and projecting obs termini to see retreat history
def read_termini(filename, year):
    """Make a dictionary of terminus positions, indexed by MEaSUREs ID. Outputs dictionary"""
    print 'Reading in MEaSUREs terminus positions for year ' + str(year)
    sf = shapefile.Reader(filename)
    fields = sf.fields[1:] #excluding the mute "DeletionFlag"
    field_names = [field[0] for field in fields]
    term_recs = sf.shapeRecords()
    termpts_dict = {}
    for r in term_recs:
        atr = dict(zip(field_names, r.record)) #dictionary of shapefile fields, so we can access GlacierID by name rather than index.  Index changes in later years.
        key = atr['GlacierID'] #MEaSUREs ID number for the glacier, found by name rather than index
        termpts_dict[key] = np.asarray(r.shape.points) #save points spanning terminus to dictionary
    return termpts_dict

# Finding intersection of terminus points with mainline--modified from hubbard-mainline-advance-v2.py
def projected_term_obs(termset, linestr):
    '''Given a termset from file input and LineString representation of a flowline, termline constructs a Shapely LineString of the terminus and returns the intersection of the two'''
    termarr = np.array(termset)
    termline = LineString(termarr)
    centrpt = termline.centroid
    arcdist = linestr.project(centrpt)
    if arcdist>0:
        return arcdist/1000
    else:
        near = linestr.distance(termline)  #in case terminus listed in MEaSUREs is farther advanced than max seaward extent of saved flowline
        return -near/1000
    
def advterm(termset, linestr):
    '''Given termset and LineString representation of a flowline, advterm finds which terminus position projects most advanced along central flowline and returns its arclength position'''
    x_term = termset[:, 0]  #Note need to change from [:, 1] to [:, 0] for x-coord, due to different data format for Hubbard
    y_term = termset[:, 1]
    projections = []
    for i in xrange(len(x_term)):
        proji = linestr.project(Point(x_term[i], y_term[i]))
        projections.append(proji)
    termmax = min(projections) #need minimum rather than max here because we are interested in the most advanced, i.e. lowest arc-length value projection of terminus
    return termmax/1000

def retterm(termset, linestr):
    '''Given termset (from file input above), retterm finds which terminus position projects most retreated (rel. 2007 terminus) along central flowline and returns its arclength position'''
    x_term = termset[:, 0]
    y_term = termset[:, 1]
    projections = []
    for i in xrange(len(x_term)):
        proji = linestr.project(Point(x_term[i], y_term[i]))
        projections.append(proji)
    termmin = max(projections) #max is the most retreated, i.e. highest arc-length value projection of terminus
    return termmin/1000

gl_termpos_fldr = 'Documents/GitHub/Data_unsynced/MEaSUREs-termini'
basefiles = ['/termini_0607_v01_2', '/termini_0708_v01_2', '/termini_0809_v01_2', '/termini_1213_v01_2', '/termini_1415_v01_2', '/termini_1516_v01_2']
obs_years = [2006, 2007, 2008, 2012, 2014, 2015] #compare with term of hindcast, 2006-2014

termini = {}
for i,b in enumerate(basefiles):
    yr = obs_years[i]
    fn = gl_termpos_fldr+b
    termini[yr] = read_termini(fn, yr) #creating dictionary for each year
    print len(termini[yr])

nw_base_fpath = 'Documents/1. Research/2. Flowline networks/Auto_selected-networks/Gld-autonetwork-GID'
projected_termini = {gid: [] for gid in glaciers_to_plot}

for gid in glaciers_to_plot:
    print 'Reading in glacier ID: '+str(gid)
    #if gid in added_jan19:
    #    filename = nw_base_fpath+str(gid)+'-date_2019-01-10.csv'
    if gid<160:
        filename = nw_base_fpath+str(gid)+'-date_2018-10-03.csv'
    else:
        filename = nw_base_fpath+str(gid)+'-date_2018-10-04.csv' #workaround because I ran these in batches and saved them with the date
    
    coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
    mainline = LineString(coords_list[0])
    for yr in obs_years:
        termpts = termini[yr][gid] #get terminus points for each year
        t = projected_term_obs(termpts, mainline) #project onto main flowline
        r = retterm(termpts, mainline) #find most retreated point
        a = advterm(termpts, mainline) #find most advanced point
        projected_termini[gid].append((a, t, r)) #add these to dictionary of projected termini per glacier

    
    

###--------------------------------------
#### PLOTTING
###--------------------------------------   
  
## Settings for plots   
labels = [str(g[1]) for g in gids_by_name if g[0] in glaciers_to_plot] #set what the glaciers will be called in plotting.  Default is simply their MEaSUREs ID
markers = ['o', '.', ',', '^', 'd', '*']
styles = ['-', ':', '-.', '-', '-', '-']
cmap = cm.get_cmap('winter')
scenario_colors = cm.get_cmap('Blues')([0.1, 0.3, 0.5, 0.7, 0.9])
#colors = cmap([0.1, 0.2, 0.3, 0.5, 0.7, 0.9])
colors = cmap(linspace(0.1, 0.9, num=len(glaciers_to_plot)))
alt_colors = cm.get_cmap('Greys')([0.2, 0.3, 0.5, 0.7, 0.9])
plot_years = 2006+np.array(testyears)


#
####terminus
plt.figure()
for j, gid in enumerate(glaciers_to_plot):
    print gid
    term_positions = full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::]
    ms_selection = mod(j, len(styles))
    plt.plot(testyears, -0.001*np.array(term_positions), linewidth=2, color='Gainsboro', linestyle=styles[ms_selection], label=labels[j])
    #plt.plot(testyears[::4], -0.001*np.array(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::])[::4], linewidth=0, marker=markers[ms_selection], ms=10, color=colors[j])
plt.legend(loc='lower left')
plt.axes().set_xlabel('Year of simulation', size=20)
plt.axes().set_ylabel('Terminus change [km]', size=20)
plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
#plt.axes().set_ylim(-100, 1)
#plt.axes().set_yticks([-75, -50, -25, 0])
plt.title('Terminus retreat of {} Greenland outlet glaciers 2006-2014 ERA-I, Tice=-10 C'.format(len(glaciers_to_plot)), fontsize=20)
plt.show()

####SINGLE NETWORK - termini vs obs
for j, gid in enumerate(glaciers_to_plot):
    sim_termini = full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'][1::]
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    plt.figure('Main line terminus change, GID{}'.format(gid))
    plt.plot(plot_years, -0.001*np.array(sim_termini), linewidth=2, color='k', linestyle='-', label='Modelled')
    plt.errorbar(obs_years, -1*obs_term_centr, yerr = e, fmt='D')
    plt.axes().set_xlabel('Year', size=30)
    plt.axes().set_ylabel('Terminus change [km]', size=30)
    plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
    plt.axes().set_xlim(2006, 2014.5)
    plt.axes().set_xticks([2006, 2008, 2010, 2012, 2014])
    plt.axes().set_xticklabels(['2006', '2008', '2010', '2012', '2014'])
    #plt.axes().set_ylim(-50, 1)
    #plt.axes().set_yticks([-50, -25, 0])
    plt.show()


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
#for j, gid in enumerate(glaciers_to_plot):
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
#plt.figure(figsize=(12,8))
#for j, gid in enumerate(glaciers_to_plot):
#    #if gid!=10:
#    ms_selection = mod(gid, len(styles))
#    plt.plot(testyears[::], scenario_sle[j], linewidth=4, color=colors[ms_selection], label=gid)
#    plt.plot(testyears[::5], scenario_sle[j][::5], linewidth=0, marker=markers[ms_selection], ms=10, color=colors[ms_selection])
#    if j==0:
#        plt.fill_between(testyears[::], y1=scenario_sle[j], y2=0, color=colors[ms_selection], alpha=0.7)  
#    else:
#        plt.fill_between(testyears[::], y1=scenario_sle[j], y2=scenario_sle[j-1], color=colors[ms_selection], alpha=0.7)     
##plt.plot([0, 20, 40, 60], [0, 14, 28, 42], color='k', linewidth=1, ls='-', alpha=0.8) #GRACE linear trend
##rect = mpatches.Rectangle((98,8.5), width=2, height=4.6, color='k', alpha=0.7) # Nick et al total projection for 2100, A1B
##rect2 = mpatches.Rectangle((98,11.3), width=2, height=6.2, color='k', alpha=0.7) # Nick et al total projection for 2100, RCP8.5 
##plt.axes().add_patch(rect)
##plt.axes().add_patch(rect2)
#plt.legend(loc='upper left')
#plt.axes().set_xlabel('Year of simulation', size=20)
#plt.axes().set_ylabel('Cumulative sea level contribution [mm]', size=20)
#plt.axes().tick_params(axis='both', length=5, width=2, labelsize=20)
##plt.axes().set_xlim(0, 100)
##plt.axes().set_xticks([0, 25, 50, 75, 100])
##plt.axes().set_ylim(0, 12)
##plt.axes().set_yticks([0, 2, 4, 6, 8, 10, 12])
#plt.show()