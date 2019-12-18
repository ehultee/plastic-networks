## Scatter plot of year-by-year dL/dt observed versus simulated
## 26 Nov 2019  EHU

import numpy as np
import matplotlib.pyplot as plt
import csv
#import collections
import glob
#from matplotlib.colors import LogNorm
from matplotlib import cm
from matplotlib.patches import Rectangle
from matplotlib.collections import PatchCollection
## Special import for SERMeQ modules
import sys
sys.path.insert(0, 'Documents/GitHub/plastic-networks')
from SERMeQ.plastic_utilities_v2 import *
from SERMeQ.GL_model_tools import *
from SERMeQ.flowline_class_hierarchy import *



###--------------------------------------
#### GLACIERS TO PLOT
###--------------------------------------
## Which glaciers are available
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
#glaciers_to_plot = [g for g in glacier_ids if g in (3, 137, 175)]
glaciers_to_plot=glacier_ids

testyears = arange(0, 9, step=0.25)#array of the years tested, with year "0" reflecting initial nominal date of MEaSUREs read-in (generally 2006)
scenarios = ('persistence',)

tempmarker = 'min10Cice' #and temperature of ice
timestepmarker = '8a_dt025a' #and total time and timestep

full_output_dicts = {}

for s in scenarios:
    scenario_output = {'Testyears': testyears}
    for gid in glaciers_to_plot:
        fn = glob.glob('Documents/GitHub/Data_unsynced/Hindcasted_networks/GID{}-*-{}-{}-{}.pickle'.format(gid, s, tempmarker, timestepmarker))[0] #using glob * to select files of multiple run dates
        lightload(fn, glacier_name = 'GID{}'.format(gid), output_dictionary = scenario_output)
    full_output_dicts[s] = scenario_output #add output from this scenario to the dictionary of all output, with scenario name as key

gl_termpos_fldr = 'Documents/GitHub/Data_unsynced/MEaSUREs-termini'
basefiles = ['/termini_0607_v01_2', '/termini_0708_v01_2', '/termini_0809_v01_2', '/termini_1213_v01_2', '/termini_1415_v01_2']
obs_years = [2006.75, 2007.75, 2009.0, 2013.0, 2014.75] #compare with years of hindcast, 2006-2014

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
    if gid in added_jan19:
        filename = nw_base_fpath+str(gid)+'-date_2019-01-10.csv'
    elif gid<160:
        filename = nw_base_fpath+str(gid)+'-date_2018-10-03.csv'
    else:
        filename = nw_base_fpath+str(gid)+'-date_2018-10-04.csv' #workaround because I ran these in batches and saved them with the date
    
    coords_list = Flowline_CSV(filename, has_width=True, flip_order=False)
    mainline = LineString(coords_list[0])
    for yr in obs_years:
        try:
            termpts = termini[yr][gid] #get terminus points for each year
            t = projected_term_obs(termpts, mainline) #project onto main flowline
            r = retterm(termpts, mainline) #find most retreated point
            a = advterm(termpts, mainline) #find most advanced point
            projected_termini[gid].append((a, t, r)) #add these to dictionary of projected termini per glacier
        except KeyError:
            print 'No terminus found in year {} for GID {}.'.format(yr, gid)
            projected_termini[gid].append((0, np.nan, 0))

###--------------------------------------
#### PLOTTING
###-------------------------------------- 

## Settings for plots   
plot_years = 2006+np.array(testyears)
ids = [i for i in range(len(plot_years)) if plot_years[i] in obs_years] #which terminus positions will we compare
yr_colors = cm.get_cmap('Blues')(linspace(0.2, 0.9, num=len(ids)))


## Compare sim vs obs terminus position
plt.figure('Yearly terminus comparison')
for j, gid in enumerate(glaciers_to_plot):
    sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    plt.errorbar(-1*obs_term_centr, -0.001*np.array(sim_termini), xerr=e, mfc='b', ecolor='b', fmt='D')
plt.plot(range(-20,2), range(-20,2), c='k', ls='-.')
plt.axes().set_aspect(1)
plt.show()


### Test from PatchCollection to make filled rectangles
def make_error_boxes(ax, xdata, ydata, xerror, yerror, colorscheme_indices,
                     edgecolor='None', barcolor='None', alpha=0.5):
    """Make a PatchCollection of filled rectangles and add it to axes.
    Keyword args same as for mpl.patches.Rectangle, except 'colorscheme_indices'
    colorscheme_indices: array of same length as xdata, ydata setting assignment of facecolors"""
    
    # Create list for all the error patches
    errorboxes = []
    # Loop over data points; create box from errors at each point
    for x, y, xe, ye, c in zip(xdata, ydata, xerror.T, yerror.T, colorscheme_indices):
        rect = Rectangle((x - xe[0], y - ye[0]), xe.sum(), ye.sum(), facecolor=cm.get_cmap('viridis_r')((c-min(colorscheme_indices))/(max(colorscheme_indices)-min(colorscheme_indices))), alpha=alpha, edgecolor=edgecolor)
        errorboxes.append(rect)
    # Create patch collection with specified colour/alpha
    pc = PatchCollection(errorboxes, match_original=True)
    # Add collection to axes
    ax.add_collection(pc)
    # Plot errorbars
    artists = ax.errorbar(xdata, ydata, xerr=xerror, yerr=yerror,
                          fmt='None', ecolor=barcolor)
    return artists


# Create figure and axes
fig, ax = plt.subplots(1)
# Call function to create error boxes
for j, gid in enumerate(glaciers_to_plot):
    sim_termini = np.take(full_output_dicts['persistence']['GID{}'.format(gid)][0]['Termini'], indices=ids)
    obs_termini = np.asarray(projected_termini[gid]) #will be of shape (len(obs_years), 3) with an entry (lower, centroid, upper) for each year
    obs_term_centr = obs_termini[:,1]
    e = np.asarray([(min(ot[0]-ot[1], ot[0]), ot[1]-ot[2]) for ot in obs_termini]).T #error lower (more advanced), upper (less advanced)
    _ = make_error_boxes(ax, -1*obs_term_centr, -0.001*np.array(sim_termini), xerror=e, yerror=0.1*np.ones(shape(e)), colorscheme_indices=obs_years)
ax.plot(range(-20,2), range(-20,2), c='k', ls='-.')
ax.set_aspect(1)
plt.show()
